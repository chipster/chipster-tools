Chipster Tools Repository
-------------------------

This repository contains the analysis tool scripts for the Chipster
software https://github.com/chipster/chipster.

For license, please see LICENSE.TXT.


Formatting
----------

R scripts are formatted with Air (https://posit-dev.github.io/air/), with
the settings in air.toml. The Format workflow on GitHub checks them with
the Air version pinned in .github/workflows/format.yml.

- VS Code: install the Air extension (Posit.air-vscode) to format on save.
- Command line: air format FILE... (or air format --check . to only check).
- Pre-commit hook: git config core.hooksPath .githooks
  The hook checks the staged R scripts, and skips the check if Air isn't
  installed.
