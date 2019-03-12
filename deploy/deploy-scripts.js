const targz = require("targz");
const fs = require("fs");
const execFile = require("child_process").execFile;
const homedir = require("os").homedir();
const path = require("path");

/**
 * Deploy tool scripts to Chipster running in OKD
 *
 * There are two modes
 *  1. Default: Deploy all tools and exit
 *  2. --watch: Watch tool script changes and deploy changed files one by one
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
      "configuration key not found: " + key + " from " + confPath
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

function main(args, confPath) {
  fs.readFile(confPath, "utf8", (err, data) => {
    if (err) {
      if (err.code == "ENOENT") {
        console.log("configuration file not found: " + confPath);
        process.exit(1);
      } else {
        throw err;
      }
    }
    try {
      var obj = JSON.parse(data);
    } catch (err) {
      console.log("configuration file parsing failed", confPath, "\n", err);
      process.exit(1);
    }

    var host = getConfiguration(obj, "okdHost");
    var token = getConfiguration(obj, "okdToken");
    var project = getConfiguration(obj, "okdProject");
    var subproject = getConfiguration(obj, "subproject");

    login(host, token, () => {
      checkProject(project, () => {
        if (args.includes("--watch")) {
          var dir = "tools";
          console.log("Watching file changes in directory '" + dir + "'...");
          for (dir of getDirectoriesRecursive(dir)) {
            watchDir(dir, filename => {
              console.log("Reload file " + filename);
              // convert windows paths
              filename = path.relative("", filename).replace(/\\/g, "/");
              console.log("Remote path " + filename);
              reloadFile(filename, subproject, () => {});
            });
          }
        } else {
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
      });
    });
  });
}

let confPath = path.join(homedir, ".chipster", "deploy-scripts.json");

main(process.argv, confPath);
