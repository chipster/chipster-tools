const fs = require("fs");
const path = require("path");
const utils = require("./deploy-utils.js");
const { bindNodeCallback, of } = require("rxjs");
const { mergeMap } = require("rxjs/operators");

/**
 * Simple launcher script for deploy-scripts.js
 *
 * This script can be used to trigger "npm install" whenever new dependencies
 * are needed. The last applied version is stored in ~/.chipster/deploy-scripts.json.
 */

// bump up this number to trigger "npm install" on all machines
const currentNodeModulesVersion = 6;

function main(argv) {
  updateNodeModulesIfNecessary()
    .pipe(mergeMap(() => deployScripts(argv)))
    .subscribe(() => {}, err => console.log(err.message ? err.message : err));
}

function updateNodeModulesIfNecessary() {
  const nodeModulesVersionKey = "nodeModulesVersion";

  return utils.getConfig().pipe(
    mergeMap(obj => {
      console.log("Check if node_modules is up-to-date");
      if (
        obj[nodeModulesVersionKey] == null ||
        obj[nodeModulesVersionKey] < currentNodeModulesVersion
      ) {
        console.log("Update node_modules");

        // it's win32 also in 64 bit Windows
        const npm = process.platform === "win32" ? "npm.cmd" : "npm"; 

        return utils.run(npm, ["install"], "deploy").pipe(
          mergeMap(() => {
            obj[nodeModulesVersionKey] = currentNodeModulesVersion;
            return bindNodeCallback(fs.writeFile)(
              utils.getConfigPath(),
              JSON.stringify(obj, null, 2),
              "utf8"
            );
          })
        );
      } else {
        return of(null);
      }
    })
  );
}

/**
 * Call the actual script
 * @param {string[]} argv full argument array of parent process
 */
function deployScripts(argv) {
  let script = "deploy-scripts.js";
  return utils.run(
    "node",
    [path.join("deploy", script), ...argv.slice(2)],
    null
  );
}

main(process.argv);
