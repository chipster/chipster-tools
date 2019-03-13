const path = require("path");
const homedir = require("os").homedir();
const fs = require("fs");
const spawn = require("child_process").spawn;
const execFile = require("child_process").execFile;

exports.getConfigPath = function() {
  return path.join(homedir, ".chipster", "deploy-scripts.json");
};

exports.getConfig = function(callback) {
  let confPath = this.getConfigPath();

  fs.readFile(confPath, "utf8", (err, data) => {
    if (err) {
      if (err.code == "ENOENT") {
        throw new Error("configuration file not found: " + confPath);
      } else {
        throw err;
      }
    }
    try {
      var obj = JSON.parse(data);
    } catch (err) {
      throw new Error("configuration file parsing failed", confPath, "\n", err);
    }

    callback(obj);
  });
};

/**
 * Run a process
 *
 * Inherit stdin, stdout and stderr from the parent process.
 * Throw an error if the process exit code is not zero.
 *
 * @param {string} cmd
 * @param {string[]} args
 * @param {string} cwd working directory path or null to inherit
 * @callback callback called after the process ends
 */
exports.run = function(cmd, args, cwd, callback) {
  const child = spawn(cmd, args, { cwd: cwd, stdio: "inherit" });
  child.on("close", code => {
    if (code != 0) {
      throw new Error(
        cmd + " " + args.join(" ") + " exited with exit code " + code
      );
    }
    callback();
  });
};

exports.runAndGetOutput = function(cmd, args, callback) {
  execFile(cmd, args, (error, stdout, stderr) => {
    if (error) {
      console.log(stdout);
      console.log(stderr);
      throw error;
    }
    callback(stdout);
  });
};
