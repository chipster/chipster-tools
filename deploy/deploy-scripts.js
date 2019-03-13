const targz = require("targz");
const fs = require("fs");
const execFile = require("child_process").execFile;
const path = require("path");
var request = require("request");
const utils = require(path.resolve(__dirname, "./deploy-utils.js"));
const WsClient = require("chipster-cli-js/lib/ws-client.js").default;
const ChipsterUtils = require("chipster-cli-js/lib/chipster-utils.js").default;

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
 * @callback callback
 */
function makePackage(packageFile, callback) {
  console.log("Packaging...");
  targz.compress(
    {
      src: "tools",
      dest: packageFile,
      tar: {
        // follow symlinks
        dereference: true,
        dmode: 0755, // needed in windows
        fmode: 0644 // needed in windows
      }
    },
    function(err) {
      if (err) {
        throw err;
      } else {
        callback();
      }
    }
  );
}

/**
 * Reload all tools from the .tar.gz package
 *
 * @param {string} packageFile path to local .tar.gz file
 * @param {string} subproject subproject name to find the correct toolbox instance
 * @callback callback
 */
function reloadAll(packageFile, subproject, callback) {
  // chmod is neeeded only if somebody has uploaded a tar without correct file modes
  var remoteScript = `
  chmod -R ugo+rwx tools
  rm -rf tools/*
  tar -xz -C tools
  `;

  reload(remoteScript, packageFile, subproject, callback);
}

/**
 * Reload a single tool script
 *
 * @param {string} file path to file to replace under the tools directory in linux format
 * @param {string} subproject subproject name to find the correct toolbox instance
 * @callback callback
 */
function reloadFile(file, subproject, callback) {
  var remoteScript = "";
  remoteScript += "rm -f " + file + "\n";
  remoteScript += "cat - > " + file + "\n";

  reload(remoteScript, file, subproject, callback);
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
 * @callback callback called when the reload has finished
 */
function reload(remoteScript, stdinFile, subproject, callback) {
  console.log("Uploading to toolbox-" + subproject + "...");

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

  var child = execFile("oc", [
    "rsh",
    "dc/toolbox-" + subproject,
    "bash",
    "-e",
    "-c",
    remoteScript
  ]);

  var reloadStarted = false;
  var logTailStarted = false;

  child.stdout.on("data", function(data) {
    var stop = false;
    // read data line by line to react to specific log messages
    // remove the last newline before splitting, because that would create an extra empty
    // line after each data block
    var lines = data.replace(/\n$/, "").split("\n");
    for (var line of lines) {
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
        process.stdout.write(line + "\n");
      }

      if (stop) {
        // stop following the log
        child.kill("SIGINT");
      }
    }
  });
  child.stderr.on("data", function(data) {
    process.stderr.write(data);
  });
  child.on("close", function(code, signal) {
    callback();
  });

  child.on("error", function(err) {
    throw err;
  });

  var stdinStream = fs.createReadStream(stdinFile);
  stdinStream.pipe(child.stdin);
}

/**
 * Watch all file changes in all directories under the directory
 *
 * Wait 100 milliseconds to see if there are more changes to come. When no
 * new changes are noticed during that time, a callback is called with the
 * filename of the last change.
 *
 * @param {string} dir
 * @callback callback
 */
function watchDir(dir, callback) {
  fs.watch(dir, (event, filename) => {
    if (filename) {
      debounce(() => {
        callback(path.resolve(dir, filename));
      }, 100);
    }
  });
}

let timeout;
/**
 * Call the function when this method hasn't been called again in the wait time
 *
 * TODO get rid of the global timeout variable
 *
 * @callback func
 * @param {number} wait in milliseconds
 */
function debounce(func, wait) {
  if (timeout) {
    // cancel previous update, because of new one
    clearTimeout(timeout);
  }
  timeout = setTimeout(() => {
    timeout = null;
    func();
  }, wait);
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
 * @callback callback
 */
function changeProject(project, callback) {
  console.log("Change to OKD project " + project);

  const child = execFile(
    "oc",
    ["project", project],
    (error, stdout, stderr) => {
      if (error) {
        console.log(stdout);
        console.log(stderr);
        throw error;
      }
      callback();
    }
  );
}

/**
 * Login to OKD if necessary
 *
 * @param {string} host
 * @param {string} token
 * @callback callback
 */
function login(host, token, callback) {
  console.log("Check if logged in to oc");

  // use "whoami" to check if logged in already (took about 1 second vs. 3 seconds to login)
  execFile("oc", ["whoami"], (error, stdout, stderr) => {
    if (error) {
      console.log("Login to " + host);
      execFile(
        "oc",
        ["login", host, "--token=" + token],
        (error, stdout, stderr) => {
          if (error) {
            console.log(stdout);
            console.log(stderr);
            throw error;
          }
          callback();
        }
      );
    }
    callback();
  });
}

/**
 * Change to a OKD project if necessary
 *
 * Store the original project in the global variable "originalProject".
 *
 * @param {string} project
 * @callback callback
 */
function checkProject(project, callback) {
  console.log("Check current project in oc");

  // check first to see if change is necessary (took 0.2 seconds vs. 2 seconds to change it)
  const child = execFile("oc", ["project", "-q"], (error, stdout, stderr) => {
    if (error) {
      console.log(stdout);
      console.log(stderr);
      throw error;
    }

    let currentProject = stdout.trim();
    if (currentProject != project) {
      originalProject = currentProject;
      console.log("current project is " + currentProject);
      changeProject(project, callback);
    } else {
      callback();
    }
  });
}

function watchAndReload(subproject, callback) {
  var dir = "tools";
  console.log("Watching file changes in directory '" + dir + "'...");
  for (dir of getDirectoriesRecursive(dir)) {
    watchDir(dir, filename => {
      console.log("Reload file " + filename);
      // convert windows paths
      filename = path.relative("", filename).replace(/\\/g, "/");
      console.log("Remote path " + filename);
      reloadFile(filename, subproject, () => callback(filename));
    });
  }
}

function packageAndReloadAll(subproject) {
  let packageDir = "build";
  if (!fs.existsSync(packageDir)) {
    fs.mkdirSync(packageDir);
  }
  let packageFile = path.join(packageDir, "tools.tar.gz");

  makePackage(packageFile, () => {
    reloadAll(packageFile, subproject, () => {
      fs.unlink(packageFile, err => {
        if (err) {
          throw err;
        }
      });
      if (originalProject != null) {
        changeProject(originalProject, () => {});
      }
    });
  });
}

function getChipsterPassword(subproject, chipsterUsername, callback) {
  console.log("Get Chipster password");
  utils.runAndGetOutput(
    "oc",
    [
      "rsh",
      "dc/auth-" + subproject,
      "bash",
      "-e",
      "-c",
      "set -o pipefail; cat security/users | grep '^" +
        chipsterUsername +
        ":' | cut -d ':' -f 2"
    ],
    password => {
      callback(password.trim());
    }
  );
}

function getChipsterToken(subproject, chipsterUsername, project, callback) {
  getChipsterPassword(subproject, chipsterUsername, password => {
    console.log("Log in to Chipster as " + chipsterUsername);
    chipsterRequest(
      "post",
      "auth",
      "tokens",
      subproject,
      project,
      chipsterUsername,
      password,
      null,
      resp => {
        callback(JSON.parse(resp).tokenKey);
      }
    );
  });
}

function chipsterRequest(
  method,
  service,
  path,
  subproject,
  project,
  username,
  password,
  body,
  callback
) {
  let uri =
    "https://" +
    service +
    "-" +
    subproject +
    "-" +
    project +
    ".rahtiapp.fi/" +
    path;
  request[method](
    uri,
    {
      auth: {
        user: username,
        pass: password
      },
      json: body != null,
      body: body
    },
    (error, response, body) => {
      if (!error && response.statusCode >= 200 && response.statusCode <= 299) {
        callback(body);
      } else {
        msg = "Chipster request failed " + method + " " + uri + "\n";
        if (response != null) {
          msg += " http status: " + response.statusCode + "\n";
        }

        if (error != null) {
          msg += " error: " + error + "\n";
        }

        if (body != null) {
          msg += " body: " + error + "\n";
        }
        throw new Error(msg);
      }
    }
  );
}

function getLatestSession(subproject, project, token, callback) {
  chipsterRequest(
    "get",
    "session-db",
    "sessions",
    subproject,
    project,
    "token",
    token,
    null,
    resp => {
      let sessions = JSON.parse(resp);

      sessions.sort(function compare(a, b) {
        var dateA = new Date(a.accessed);
        var dateB = new Date(b.accessed);
        return dateB - dateA;
      });
      callback(sessions[0]);
    }
  );
}

function getLatestJob(subproject, project, token, sessionId, callback) {
  chipsterRequest(
    "get",
    "session-db",
    "sessions/" + sessionId + "/jobs",
    subproject,
    project,
    "token",
    token,
    null,
    resp => {
      let jobs = JSON.parse(resp);

      jobs.sort(function compare(a, b) {
        var dateA = new Date(a.created);
        var dateB = new Date(b.created);
        return dateB - dateA;
      });

      callback(jobs[0]);
    }
  );
}

function runJob(subproject, project, token, sessionId, job, callback) {
  chipsterRequest(
    "post",
    "session-db",
    "sessions/" + sessionId + "/jobs",
    subproject,
    project,
    "token",
    token,
    job,
    resp => {
      console.log("job created", resp);
      callback(resp.jobId);
    }
  );
}

function followJob(sessionId, jobId, token, subproject, job) {
  ChipsterUtils.getRestClient(
    "https://" + subproject + "-" + project + ".rahtiapp.fi",
    token
  ).subscribe(
    restClient => {
      let wsClient = new WsClient(restClient);
      wsClient.connect(sessionId);

      wsClient.getJobState$(jobId).subscribe(
        job => {
          console.log("*", job.state, "(" + (job.stateDetail || "") + ")");
        },
        err => console.error("failed to get the job state", err),
        () => {
          wsClient.disconnect();
        }
      );

      wsClient.getJobScreenOutput$(jobId).subscribe(
        output => {
          process.stdout.write(output);
        },
        err => console.debug("job screen output error", err)
      );
      wsClient.getJobOutputDatasets$(jobId).subscribe(
        dataset => {
          console.log(
            "* dataset created: " + dataset.name.padEnd(24) + dataset.datasetId
          );
        },
        err => console.debug("job output file error", err)
      );
    },
    err => console.log("rest client error", err)
  );
}

function watchAndRun(subproject, chipsterUserId, project, callback) {
  let chipsterUsername = chipsterUserId.split("/")[1];
  getChipsterToken(subproject, chipsterUsername, project, token => {
    watchAndReload(subproject, () => {
      getLatestSession(subproject, project, token, session => {
        getLatestJob(subproject, project, token, session.sessionId, job => {
          job.state = "NEW";
          job.jobId = null;
          console.log(
            "Rerun latest job " + job.toolId + " in session " + session.name
          );
          runJob(subproject, project, token, session.sessionId, job, jobId => {
            followJob(session.sessionId, jobId, token, subproject, project);
          });
        });
      });
    });
  });
}

function main(args) {
  utils.getConfig(obj => {
    var host = getConfiguration(obj, "okdHost");
    var token = getConfiguration(obj, "okdToken");
    var project = getConfiguration(obj, "okdProject");
    var subproject = getConfiguration(obj, "subproject");
    var chipsterUserId = getConfiguration(obj, "chipsterUsername");

    login(host, token, () => {
      checkProject(project, () => {
        if (args.includes("--watch")) {
          watchAndReload(subproject, () => {});
        } else if (args.includes("--run")) {
          watchAndRun(subproject, chipsterUserId, project, () => {});
        } else {
          packageAndReloadAll(subproject);
        }
      });
    });
  });
}

main(process.argv);
