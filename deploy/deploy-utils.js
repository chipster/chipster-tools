const path = require("path");
const homedir = require("os").homedir();
const fs = require("fs");
const spawn = require("child_process").spawn;
const execFile = require("child_process").execFile;
const { Subject, of, bindCallback } = require("rxjs");
const { map, mergeMap } = require("rxjs/operators");

exports.getConfigPath = function() {
  return path.join(homedir, ".chipster", "deploy-scripts.json");
};

exports.getConfig = function() {
  let confPath = this.getConfigPath();
  let subject = new Subject();

  fs.readFile(confPath, "utf8", (err, data) => {
    if (err) {
      if (err.code == "ENOENT") {
        subject.error(new Error("configuration file not found: " + confPath));
      } else {
        subject.error(err);
      }
    } else {
      try {
        let parsed = JSON.parse(data);
        subject.next(parsed);
        subject.complete();
      } catch (err) {
        console.log("conf parse error", err);
        subject.error(
          new Error("configuration file parsing failed", confPath, "\n", err)
        );
      }
    }
  });
  return subject;
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
 */
exports.run = function(cmd, args, cwd) {
  return of(null).pipe(
    mergeMap(() => {
      // how would we convert this to an observable with bindCallback, when we are
      // are interested only about the "close" event?
      let subject = new Subject();
      const child = spawn(cmd, args, { cwd: cwd, stdio: "inherit" });
      child.on("close", code => {
        if (code != 0) {
          subject.error(
            new Error(
              cmd + " " + args.join(" ") + " exited with exit code " + code
            )
          );
        }
        subject.next();
        subject.complete();
      });
      return subject;
    })
  );
};

exports.runAndGetOutput = function(cmd, args) {
  // use regular bindCallback (instead of bindNodeCallback) to access stdout
  // and stderr also when there is an error
  return bindCallback(execFile)(cmd, args).pipe(
    map(res => {
      let err = res[0];
      let stdout = res[1];
      let stderr = res[2];

      if (err) {
        console.log(stdout);
        console.log(stderr);
        throw err;
      }
      return stdout;
    })
  );
};
