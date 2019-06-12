const targz = require("targz");
const fs = require("fs");
const execFile = require("child_process").execFile;
const path = require("path");
var request = require("request");
const utils = require(path.resolve(__dirname, "./deploy-utils.js"));
const WsClient = require("chipster-cli-js/lib/ws-client.js").default;
const RestClient = require("chipster-nodejs-core/src/rest-client.js")
  .RestClient;
const { of, bindNodeCallback, Subject, merge, forkJoin } = require("rxjs");
const {
  mergeMap,
  catchError,
  debounceTime,
  map,
  tap,
  toArray,
  takeUntil,
  take
} = require("rxjs/operators");

/**
 * Deploy tool scripts to Chipster running in OKD
 *
 * There are different modes
 *  1. Default: Deploy all tools and exit
 *  2. --watch: Watch tool script changes and deploy changed files one by one
 *  3. --run: Watch tool script changes and deploy changed files and restart the latest job
 */

let originalProject;

function flatten(lists) {
  return lists.reduce((a, b) => a.concat(b), []);
}

function getDirectories(srcpath) {
  return fs
    .readdirSync(srcpath)
    .map(file => path.join(srcpath, file))
    .filter(path => fs.statSync(path).isDirectory());
}

/**
 * Get all directories under the srcpath
 *
 * https://stackoverflow.com/a/40896897
 *
 * @param {string} srcpath
 */
function getDirectoriesRecursive(srcpath) {
  return [
    srcpath,
    ...flatten(getDirectories(srcpath).map(getDirectoriesRecursive))
  ];
}

/**
 * Package all files in the directory "tools" to a .tar.gz file
 *
 * @param {string} packageFile
 */
function makePackage(packageFile) {
  console.log("Packaging...");
  return bindNodeCallback(targz.compress)({
    src: "tools",
    dest: packageFile,
    tar: {
      // follow symlinks
      dereference: true,
      dmode: 0755, // needed in windows
      fmode: 0644 // needed in windows
    }
  });
}

/**
 * Reload all tools from the .tar.gz package
 *
 * @param {string} packageFile path to local .tar.gz file
 * @param {string} subproject subproject name to find the correct toolbox instance
 */
function reloadAll(packageFile, subproject) {
  // chmod is neeeded only if somebody has uploaded a tar without correct file modes
  var remoteScript = `
  rm -rf tools
  mkdir tools
  tar -xz -C tools
  `;

  return reload(remoteScript, packageFile, subproject);
}

/**
 * Reload a single tool script
 *
 * @param {string} file path to file to replace under the tools directory in linux format
 * @param {string} subproject subproject name to find the correct toolbox instance
 */
function reloadFile(file, subproject, onlyErrors) {
  var remoteScript = "";
  remoteScript += "rm -f " + file + "\n";
  remoteScript += "cat - > " + file + "\n";

  return reload(remoteScript, file, subproject, onlyErrors);
}

/**
 * Reload tools in the toolbox of the given subproject
 *
 * Tool scripts can be updated using the stdinFile and remoteScript before the
 * reload.
 *
 * After the reload is triggered, we try to parse the relevant section from the
 * toolbox log.
 *
 * @param {string} remoteScript bash script executed in the toolbox container
 * @param {string} stdinFile path to local file that is passed to the stdin of the remoteScript
 * @param {string} subproject name to find the correct toolbox instance
 */
function reload(remoteScript, stdinFile, subproject, onlyErrors = false) {
  console.log("Reload toolbox-" + subproject);

  logTailStartsMark = "LOG TAIL STARTS";

  remoteScript +=
    `
  echo "Reloading..."
  touch .reload/touch-me-to-reload-tools
  echo ` +
    logTailStartsMark +
    `
  tail -f logs/chipster.log & sleep 5; kill %%
  echo Done
  `;

  let child = execFile("oc", [
    "rsh",
    "dc/toolbox-" + subproject,
    "bash",
    "-e",
    "-c",
    remoteScript
  ]);

  let reloadStarted = false;
  let logTailStarted = false;
  let onlyErrorsBuffer = "";

  let subject = new Subject();

  child.stdout.on("data", data => {
    let stop = false;
    // read data line by line to react to specific log messages
    // remove the last newline before splitting, because that would create an extra empty
    // line after each data block
    let lines = data.replace(/\n$/, "").split("\n");
    for (let line of lines) {
      if (line.indexOf(logTailStartsMark) != -1) {
        logTailStarted = true;
      }

      // skip messages until the reload is started
      if (line.indexOf("tool reload requested") != -1) {
        reloadStarted = true;
      }
      if (reloadStarted && line.indexOf("tools reload done") != -1) {
        stop = true;
      }

      if (reloadStarted && line.indexOf("failed to reload tools") != -1) {
        stop = true;
      }

      // show output until log tailing starts. Then skip all old log messages until the tool reload starts
      if (!logTailStarted || reloadStarted) {
        if (onlyErrors) {
          onlyErrorsBuffer += line + "\n";
        } else {
          process.stdout.write(line + "\n");
        }
      }

      if (stop) {
        // stop following the log
        child.kill("SIGINT");
      }
    }
  });
  child.stderr.on("data", data => {
    // silence these repeating messages
    if (
      data.indexOf("Defaulting container name to ") == -1 &&
      data.indexOf(" to see all of the containers in this pod.") == -1
    ) {
      process.stderr.write(data);
    }
  });
  child.on("close", (code, signal) => {
    subject.next();
    subject.complete();
  });

  child.on("error", err => {
    if (onlyErrors) {
      process.stdout.write(onlyErrorsBuffer);
    }
  });

  var stdinStream = fs.createReadStream(stdinFile);
  stdinStream.pipe(child.stdin);

  return subject;
}

/**
 * Watch all file changes in all directories under the directory
 *
 * @param {string} dir
 */
function watchDir(dir) {
  let subject = new Subject();
  fs.watch(dir, (event, filename) => {
    if (filename) {
      let absolutePath = path.resolve(dir, filename);
      let relativePath = path.relative("", absolutePath);
      subject.next(relativePath);
    }
  });
  return subject;
}

/**
 * Get a value of the field in object
 *
 * Throw error if the field is not found (or the value is null).
 *
 * @param {Object} obj
 * @param {string} key
 */
function getConfiguration(obj, key) {
  if (obj[key] == null) {
    throw new Error(
      "configuration key not found: " + key + " from " + utils.getConfigPath()
    );
  }
  return obj[key];
}

/**
 * Change to given OKD project
 *
 * @param {string} project
 */
function changeProject(project) {
  console.log("Change to OKD project " + project);
  return utils.run("oc", ["project", project]);
}

/**
 * Login to OKD if necessary
 *
 * @param {string} host
 * @param {string} token
 */
function login(host, token) {
  console.log("Check if logged in to oc");

  // use "whoami" to check if logged in already (took about 1 second vs. 3 seconds to login)
  return utils.runAndGetOutput("oc", ["whoami"]).pipe(
    catchError(err => {
      console.log("Log in");
      return utils.run("oc", ["login", host, "--token=" + token]);
    })
  );
}

/**
 * Change to a OKD project if necessary
 *
 * Store the original project in the global variable "originalProject".
 *
 * @param {string} project
 */
function checkProject(project) {
  console.log("Check current project in oc");

  // check first to see if change is necessary (took 0.2 seconds vs. 2 seconds to change it)
  return utils.runAndGetOutput("oc", ["project", "-q"]).pipe(
    mergeMap(stdout => {
      let currentProject = stdout.trim();
      if (currentProject != project) {
        originalProject = currentProject;
        console.log("current project is " + currentProject);
        return changeProject(project);
      } else {
        return of(null);
      }
    })
  );
}

function watchAndReload(subproject, toolsDir, onlyErrors = false) {
  printWatch(toolsDir);
  let dirWatches$ = getDirectoriesRecursive(toolsDir).map(watchDir);
  return merge(...dirWatches$).pipe(
    // wait 100 ms and discard this event if there is a new one
    debounceTime(100),
    mergeMap(filename => {
      console.log("Update tool " + filename);
      // convert windows paths
      filename = filename.replace(/\\/g, "/");
      return reloadFile(filename, subproject, onlyErrors);
    })
  );
}

function packageAndReloadAll(subproject) {
  let packageDir = "build";
  let packageFile = path.join(packageDir, "tools.tar.gz");

  return bindNodeCallback(fs.access)(packageDir).pipe(
    catchError(err => {
      if (err.code == "ENOENT") {
        return bindNodeCallback(fs.mkdir)(packageDir);
      } else {
        return throwError(err);
      }
    }),
    mergeMap(() => makePackage(packageFile)),
    mergeMap(() => reloadAll(packageFile, subproject)),
    mergeMap(() => bindNodeCallback(fs.unlink)(packageFile)),
    mergeMap(() => {
      if (originalProject != null) {
        return changeProject(originalProject, () => {});
      } else {
        return of(null);
      }
    })
  );
}

function getChipsterPassword(subproject, chipsterUsername) {
  console.log("Get Chipster password");
  args = [
    "rsh",
    "dc/auth-" + subproject,
    "bash",
    "-e",
    "-c",
    "set -o pipefail; cat security/users | grep '^" +
      chipsterUsername +
      ":' | cut -d ':' -f 2"
  ];

  return utils.runAndGetOutput("oc", args).pipe(map(stdout => stdout.trim()));
}

function getLatestSession(restClient, webServerUri, chipsterUserId) {
  return restClient.getSessions().pipe(
    map(sessions => {
      const last = sortByDate(sessions, "accessed")[0];
      if (last) {
        return last;
      }
      throw new Error(
        "last session not found. Please login to " +
          webServerUri +
          " as " +
          chipsterUserId +
          " and open a session"
      );
    })
  );
}

function getLatestJob(restClient, session, webServerUri, chipsterUserId) {
  return restClient.getJobs(session.sessionId).pipe(
    map(jobs => {
      const last = sortByDate(jobs, "created")[0];
      if (last) {
        return last;
      }
      throw new Error(
        "last job not found from session '" +
          session.name +
          "' in " +
          webServerUri
      );
    })
  );
}

function sortByDate(array, field) {
  array.sort(function compare(a, b) {
    var dateA = new Date(a[field]);
    var dateB = new Date(b[field]);
    return dateB - dateA;
  });
  return array;
}

function followJob(jobId, wsClient) {
  let datasets = [];

  let screenOutput$ = wsClient
    .getJobScreenOutput$(jobId)
    // wait for stdout flush
    .pipe(tap(output => process.stdout.write(output)));

  let outputDatasets$ = wsClient.getJobOutputDatasets$(jobId).pipe(
    tap(dataset => {
      console.log(
        "* dataset created: " + dataset.name.padEnd(24) + dataset.datasetId
      );
      datasets.push(dataset);
    }),
    //TODO find out why this doesn't complete when there are no outputs, e.g. when the job fails
    takeUntil(wsClient.getJobState$(jobId).pipe(toArray()))
  );

  let jobState$ = wsClient.getJobState$(jobId).pipe(
    tap(job => {
      console.log("*", job.state, "(" + (job.stateDetail || "") + ")");
    })
  );

  return forkJoin(
    screenOutput$.pipe(
      // emit even if empty
      toArray(),
      mergeMap(() => of(null))
    ),
    outputDatasets$.pipe(
      toArray(),
      mergeMap(() => of(null))
    ),
    jobState$.pipe(
      toArray(),
      mergeMap(() => of(null))
    )
  ).pipe(
    tap(() => wsClient.disconnect()),
    map(() => {
      let jobObjects = {
        jobId: jobId,
        datasets: datasets
      };
      return jobObjects;
    })
  );
}

function getRestClientPatched(isClient, token, serviceLocatorUri) {
  let restClient = new RestClient(isClient, token, serviceLocatorUri);

  // patch this methods in RestClient, for some reason authenticated reqeusts get stuck in the default implementation
  restClient.get = (uri, headers) => {
    return httpRequest("get", uri, headers, restClient);
  };

  restClient.delete = (uri, headers) => {
    return httpRequest("delete", uri, headers, restClient);
  };

  return restClient;
}

function httpRequest(method, uri, headers, restClient) {
  let options = {
    method: method,
    uri: uri,
    headers: headers
  };
  //TODO why authenticated requests get stuck if we do this with bindNodeCallback()?
  let subject = new Subject();
  request(options, (error, response, body) => {
    if (error) {
      console.log("request error", options.method, options.uri, err);
      subject.error(err);
    }
    subject.next({
      response: response,
      body: body
    });
  });

  return subject.pipe(map(restClient.handleResponse));
}

/**
 * Like ChipsterUtils.login(), but using the patched RestClient
 *
 * @param {*} webServerUri
 * @param {*} username
 * @param {*} password
 */
function loginChipster(webServerUri, username, password) {
  // get the service locator address
  return getRestClientPatched(webServerUri)
    .getServiceLocator(webServerUri)
    .pipe(
      // get token
      mergeMap(serviceLocatorUrl => {
        let guestClient = getRestClientPatched(true, null, serviceLocatorUrl);
        return guestClient.getToken(username, password);
      })
    );
}

/**
 * * Like ChipsterUtils.getRestClient(), but using the patched RestClient
 *
 * @param {*} webServerUri
 * @param {*} chipsterToken
 */
function getRestClient(webServerUri, chipsterToken) {
  return getRestClientPatched(true, null, null)
    .getServiceLocator(webServerUri)
    .pipe(
      map(serviceLocatorUri => {
        return getRestClientPatched(true, chipsterToken, serviceLocatorUri);
      })
    );
}

function deleteJob(restClient, jobObjects) {
  if (jobObjects) {
    let msg = "Delete the previous job ";
    if (jobObjects.datasets.length > 0) {
      msg += "and its " + jobObjects.datasets.length + " output file(s)";
    }
    console.log(msg);
    requests = [];

    // complete the observables with take(1) so that forkJoin emits
    // httpRequest() doesn't complete the observable, because that is what RestClient seems to assume
    requests.push(
      restClient.deleteJob(jobObjects.sessionId, jobObjects.jobId).pipe(take(1))
    );

    jobObjects.datasets.forEach(d => {
      requests.push(
        restClient
          .deleteDataset(jobObjects.sessionId, d.datasetId)
          .pipe(take(1))
      );
    });

    return forkJoin(requests);
  } else {
    return of(null);
  }
}

function getWsClient(restClient, sessionId) {
  let wsClient = new WsClient(restClient);
  // this doesn't wait, hopefully wsClient handles this internally
  wsClient.connect(sessionId, true);
  return wsClient;
}

function runLatestJob(session, restClient, webServerUri, chipsterUserId) {
  let wsClient;
  let job2;

  return getLatestJob(restClient, session, webServerUri, chipsterUserId).pipe(
    tap(job => {
      job.state = "NEW";
      job.jobId = null;
      console.log(
        "Restart latest job " + job.toolId + " in session " + session.name
      );
      job2 = job;
    }),
    map(() => getWsClient(restClient, session.sessionId)),
    tap(wc => (wsClient = wc)),
    mergeMap(() => restClient.postJob(session.sessionId, job2)),
    mergeMap(jobId => followJob(jobId, wsClient)),
    map(o => {
      o.sessionId = session.sessionId;
      return o;
    })
  );
}

function watchAndRun(subproject, chipsterUserId, project, toolsDir) {
  let chipsterUsername = chipsterUserId.split("/")[1];
  let webServerUri = "https://" + subproject + "-" + project + ".rahtiapp.fi";
  let previousJobObjects;
  return getChipsterPassword(subproject, chipsterUsername).pipe(
    mergeMap(password => {
      return loginChipster(webServerUri, chipsterUsername, password);
    }),
    map(token => token.tokenKey),
    mergeMap(chipsterToken => {
      return getRestClient(webServerUri, chipsterToken);
    }),
    mergeMap(restClient => {
      return watchAndReload(subproject, toolsDir, true).pipe(
        mergeMap(() => deleteJob(restClient, previousJobObjects)),
        mergeMap(() =>
          getLatestSession(restClient, webServerUri, chipsterUserId)
        ),
        mergeMap(session =>
          runLatestJob(session, restClient, webServerUri, chipsterUserId)
        ),
        tap(jo => {
          previousJobObjects = jo;
        })
      );
    })
  );
}

function printWatch(dir) {
  console.log(
    "Watch file changes in directory '" + dir + "'... (Ctrl+C to interrupt)"
  );
}

function main(args) {
  let token;
  let project;
  let subproject;
  let chipsterUserId;

  toolsDir = "tools";

  utils
    .getConfig()
    .pipe(
      mergeMap(obj => {
        let host = getConfiguration(obj, "okdHost");
        token = getConfiguration(obj, "okdToken");
        project = getConfiguration(obj, "okdProject");
        subproject = getConfiguration(obj, "subproject");
        chipsterUserId = getConfiguration(obj, "chipsterUsername");
        return login(host, token);
      }),
      mergeMap(() => checkProject(project)),
      mergeMap(() => {
        if (args.includes("--watch")) {
          return watchAndReload(subproject, toolsDir).pipe(
            tap(() => printWatch(toolsDir))
          );
        } else if (args.includes("--run")) {
          return watchAndRun(
            subproject,
            chipsterUserId,
            project,
            toolsDir
          ).pipe(tap(() => printWatch(toolsDir)));
        } else {
          return packageAndReloadAll(subproject);
        }
      })
    )
    .subscribe(null, err => console.log(err));
}

main(process.argv);
