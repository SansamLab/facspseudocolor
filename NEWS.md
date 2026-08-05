# facspseudocolor 0.1.0.9000

## Quality assurance

- Corrected the Python-boundary test so source-checkout testing verifies the package-build exclusion rule while installed-package testing directly verifies that repository-only Python tools are absent.

# facspseudocolor 0.1.0

## First lab-ready release

`facspseudocolor` is now an installable R package for validating, normalizing, quantifying, and plotting event-level signal-versus-DNA flow-cytometry data. It supports EdU incorporation, protein-of-interest measurements, and phospho-histone H3 analysis.

## Analysis modes

- Added EdU analysis with an independently fitted background model for every acquisition. EdU-negative events determine the signal baseline without requiring a separate background-control sample.
- Added protein-of-interest analysis using matched background-control samples and DNA-dependent background correction.
- Added pH3 analysis using user-defined FlowJo positivity gates, per-sample G1 DNA normalization, explicit cell-cycle boundaries, unassigned-event reporting, and G2/M boundary-sensitivity diagnostics.
- Added exact FlowJo gate-geometry extraction and gate-overlay plots for pH3 analyses.
- Added whole-population and cell-cycle-phase quantitation for background-subtracted and normalized signal representations.
- Added optional within-replicate reference normalization after background subtraction.

## Plotting and reports

- Added publication-style pseudocolor plots with consistent axes across samples and replicates.
- Made background-subtracted signal the default pseudocolor and quantitation representation.
- Added automatic and manually configurable display offsets for background-subtracted pseudocolor plots. Display offsets do not affect quantitation.
- Corrected POI cutoff-line placement when a background-subtracted display offset is used.
- Added complete, pseudocolor, quantitation, cell-cycle, diagnostic, and pH3 Quarto report templates.
- Added validated presentation overrides, partial condition-color mappings, and colorblind-friendly palette options.
- Made `cowplot` the default layout. `plotgardener` remains available as an optional dependency for fixed-size publication layouts and editable figure bundles.

## Configuration and reproducibility

- Added strict YAML configuration and input validation with explicit provenance.
- Added a Shiny-enabled Quarto configurator for assigning files, samples, conditions, biological and technical replicates, analysis roles, channels, and output settings.
- Added optional Excel-based configuration for larger experiments.
- Added structured analysis, quantitation, plotting, and safe-saving APIs.
- Added durable analysis and editable figure-bundle RDS artifacts.
- Added pure in-memory EdU, POI, and pH3 normalization APIs.
- Refused accidental overwriting of analysis results unless explicitly requested.

## Installation and portability

- Converted the workflow into an installable R package requiring R 4.2 or newer.
- Made Python optional and limited it to rebuilding population exports from FlowJo files. Python is not required to install, load, test, or use the R package with existing CSV exports.
- Made `plotgardener` optional so the standard package installation does not require Bioconductor.
- Added direct installation instructions for lab users through `pak`.
- Added automated package checks on Linux, macOS, and Windows using current and previous R releases.

## Quality assurance

- Added unit, integration, numerical-regression, artifact, plotting, and Python-boundary tests.
- Added cross-platform artifact-path tests.
- Added automated GitHub pull-request checks.
- Added contributing and security policies.
