const targz = require("targz");
const fs = require("fs");
const execFile = require("child_process").execFile;
const homedir = require("os").homedir();
const path = require("path");

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

// https://stackoverflow.com/a/40896897
function getDirectoriesRecursive(srcpath) {
  return [
    srcpath,
    ...flatten(getDirectories(srcpath).map(getDirectoriesRecursive))
  ];
}

function makePackage(packageFile, callback) {
  console.log("Packaging...");
  targz.compress(
    {
      src: "tools",
      dest: packageFile,
      tar: {
        // follow symlinks
        dereference: true
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

function reloadAll(stdinFile, subproject, callback) {
  var remoteScript = `
  rm -rf tools
  mkdir tools
  pushd tools
  tar -xz
  popd
  `;

  reload(remoteScript, stdinFile, subproject, callback);
}

function reloadFile(file, subproject, callback) {
  var remoteScript = "";
  remoteScript += "rm -f " + file + "\n";
  remoteScript += "cat - > " + file + "\n";

  reload(remoteScript, file, subproject, callback);
}

function reload(remoteScript, stdinFile, subproject, callback) {
  console.log("Uploading to toolbox-" + subproject + "...");

  remoteScript += `
  echo "Reloading..."
  touch .reload/touch-me-to-reload-tools
  tail -f logs/chipster.log & sleep 5; kill %%
  echo Done
  `;

  var child = execFile("oc", [
    "rsh",
    "dc/toolbox-" + subproject,
    "bash",
    "-c",
    remoteScript
  ]);

  var reloadStarted = false;

  child.stdout.on("data", function(data) {
    var stop = false;
    // read data line by line to react to specific log messages
    // remove the last newline before splitting, because that would create an extra empty
    // line after each data block
    var lines = data.replace(/\n$/, "").split("\n");
    for (var line of lines) {
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

      if (reloadStarted) {
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

function watchDir(dir, callback) {
  fs.watch(dir, (event, filename) => {
    if (filename) {
      debounce(() => {
        callback(path.resolve(dir, filename));
      }, 100);
    }
  });
}

var timeout;
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

function getConfiguration(obj, key) {
  if (obj[key] == null) {
    throw new Error(
      "configuration key not found: " + key + " from " + confPath
    );
  }
  return obj[key];
}

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

function login(host, token, callback) {
  console.log("Check if logged in to oc");

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

function checkProject(project, callback) {
  console.log("Check current project in oc");

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
              reloadFile(filename, subproject, () => {});
            });
          }
        } else {
          let packageFile = path.join("build", "tools.tar.gz");
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
