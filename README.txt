Chipster Tools Repository
-------------------------

This repository contains the analysis tool scripts for the Chipster
software https://github.com/chipster/chipster.

For license, please see LICENSE.TXT.


Formatting
----------

R scripts are formatted with Air (https://posit-dev.github.io/air/), with
the settings in air.toml, and Python scripts with Ruff
(https://docs.astral.sh/ruff/), with the settings in ruff.toml. The Format
workflow on GitHub checks them with the versions pinned in
.github/workflows/format.yml.

- VS Code: install the Air (Posit.air-vscode) and Ruff (charliermarsh.ruff)
  extensions to format on save.
- Command line: air format FILE... and ruff format FILE... (or
  air format --check . and ruff format --check . to only check).
- Pre-commit hook: git config core.hooksPath .githooks
  The hook checks the staged R and Python scripts, and skips a check if its
  formatter isn't installed.
