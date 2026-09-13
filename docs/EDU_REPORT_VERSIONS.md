# EdU report versions

EdU report templates are retained as a series. A new design receives a new
filename instead of replacing an established template.

## Compatibility policy

- Existing report filenames remain installed so old render commands continue
  to resolve.
- A versioned report's reader-facing sections and parameter contract are frozen.
  Correctness or security fixes may be applied, but an intentional layout,
  feature, or parameter-contract change requires a new report version.
- Every supported template must continue accepting exactly one of `config` or
  `analysis_rds`.
- Package data objects retain documented compatibility aliases for their stated
  support window. A report must not silently redirect an old field name to a
  different scientific meaning.
- `inst/quarto/report-catalog.yml` is the machine-readable inventory and status
  record.

## Current series

| Report ID | Template | Status | Purpose |
|---|---|---|---|
| `edu_legacy_v1` | `facs_edu_pseudocolor_output_contract.qmd` | Stable compatibility | The established report and original command target from the 0.1 development line. |
| `edu_standard_v2` | `facs_edu_standard_v2.qmd` | Current development | The enhanced report with interactive apex lines, phase boxes, manuscript-scale PDF downloads, gating cards, phase cards, and background-correction diagnostics. |
| `edu_model_gated_experimental` | `facs_edu_model_gated.qmd` | Experimental | Compact report for the explicitly experimental frozen model-derived gates. |

The unfinished configurator is cataloged separately and is not a stable report.

## Rendering a specific version

From the package source directory:

```bash
# Established legacy-compatible report
quarto render inst/quarto/facs_edu_pseudocolor_output_contract.qmd \
  -P config=/absolute/path/to/config.yml

# Enhanced standard v2 report
quarto render inst/quarto/facs_edu_standard_v2.qmd \
  -P config=/absolute/path/to/config.yml
```

An installed package exposes the same filenames under its `quarto` directory.
Record the package version and selected report ID with each retained analysis.

## Starting the next report

Copy the current version to the next permanent filename, update the catalog,
add a section/parameter contract test for the new file, and leave all older
templates and their tests in place.
