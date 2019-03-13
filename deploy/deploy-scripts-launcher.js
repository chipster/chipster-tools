const fs = require("fs");
const path = require("path");
const utils = require("./deploy-utils.js");

/**
 * Simple launcher script for deploy-scripts.js
 *
 * This script can be used to trigger "npm install" whenever new dependencies
 * are needed. The last applied version is stored in ~/.chipster/deploy-scripts.json.
 */

// bump up this number to trigger "npm install" on all machines
const currentNodeModulesVersion = 3;

function main(argv, confPath) {
  const nodeModulesVersionKey = "nodeModulesVersion";

  utils.getConfig(obj => {
    console.log("Check if node_modules is up-to-date");
    if (
      obj[nodeModulesVersionKey] == null ||
      obj[nodeModulesVersionKey] < currentNodeModulesVersion
    ) {
      console.log("Update node_modules");
      utils.run("npm", ["install"], "deploy", () => {
        obj[nodeModulesVersionKey] = currentNodeModulesVersion;
        fs.writeFileSync(
          utils.getConfigPath(),
          JSON.stringify(obj, null, 2),
          "utf8"
        );
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
  utils.run(
    "node",
    [path.join("deploy", "deploy-scripts.js"), ...argv.slice(2)],
    null,
    () => {}
  );
}

main(process.argv);
