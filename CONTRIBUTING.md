# Contributing to facspseudocolor

## Development setup

Use R 4.2 or newer. From the repository root, install the development and
runtime dependencies, then install the package:

```r
install.packages(c(
  "cowplot", "ggplot2", "MASS", "openxlsx", "readxl", "scales",
  "shiny", "testthat", "withr", "xml2", "yaml"
))
install.packages("devtools")
devtools::install(dependencies = FALSE)
```

Python is optional and is used only to rebuild CSV exports from FlowJo files.
See `docs/PYTHON_INTERFACE.md` before changing that integration boundary.

## Before opening a pull request

Run these commands from the repository root:

```bash
Rscript -e 'devtools::document()'
Rscript -e 'devtools::test()'
R CMD build .
R CMD check --no-manual --as-cran facspseudocolor_*.tar.gz
```

Do not commit generated analysis results, local R libraries, virtual
environments, or rendered Quarto output. Add or update tests for behavioral
changes, and update `NEWS.md` for user-visible changes.

## Scientific behavior changes

Changes to normalization, gating, density calculations, or summary statistics
must include a numerical regression test or a small reproducible fixture.
Describe the scientific rationale and expected effect in the pull request.

## Pull requests

Keep each pull request focused. Explain the user-facing problem, the approach,
the validation performed, and any compatibility or scientific interpretation
concerns. All automated package checks must pass before merging.
