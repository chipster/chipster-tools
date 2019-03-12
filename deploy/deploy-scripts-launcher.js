const fs = require("fs");
const spawn = require("child_process").spawn;
const path = require("path");
const homedir = require("os").homedir();

/**
 * Simple launcher script for deploy-scripts.js
 *
 * This script can be used to trigger "npm install" whenever new dependencies
 * are needed. The last applied version is stored in ~/.chipster/deploy-scripts.json.
 */

// bump up this number to trigger "npm install" on all machines
const currentNodeModulesVersion = 1;

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
function run(cmd, args, cwd, callback) {
  const child = spawn(cmd, args, { cwd: cwd, stdio: "inherit" });
  child.on("close", code => {
    if (code != 0) {
      throw new Error(cmd + args.join(" ") + " exited with exit code " + code);
    }
    callback();
  });
}

function main(argv, confPath) {
  const nodeModulesVersionKey = "nodeModulesVersion";
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

    console.log("Check if node_moduels is up-to-date");
    if (
      obj[nodeModulesVersionKey] == null ||
      obj[nodeModulesVersionKey] < currentNodeModulesVersion
    ) {
      console.log("Update node_modules");
      run("npm", ["install"], "deploy", () => {
        obj[nodeModulesVersionKey] = currentNodeModulesVersion;
        fs.writeFileSync(confPath, JSON.stringify(obj), "utf8");
        deployScripts(argv);
      });
    } else {
      deployScripts(argv);
    }
  });
}

/**
 * Call the actual script
 * @param {string[]} argv full argument array of parent process
 */
function deployScripts(argv) {
  run(
    "node",
    [path.join("deploy", "deploy-scripts.js"), ...argv.slice(2)],
    null,
    () => {}
  );
}

let confPath = path.join(homedir, ".chipster", "deploy-scripts.json");

main(process.argv, confPath);
