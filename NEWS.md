# facspseudocolor 0.0.0.9000

- Added background-subtracted pseudocolor display with an automatic offset
  rounded from the pooled raw G1 target level to the nearest power of ten, plus
  an explicit manual-offset option. Quantitation remains unshifted.
- Added PH3 mode with user-defined FlowJo positivity, G1-based per-sample DNA
  normalization, explicit phase boundaries, all-Single-Cells denominators,
  unassigned-event reporting, and G2/M boundary-sensitivity diagnostics.
- Added a separate exact FlowJo gate-geometry extractor and a focused pH3
  Quarto report with editable gate-overlaid pseudocolor panels.
- Converted reusable EdU and POI analysis into an installable R package.
- Added strict YAML and input validation with explicit provenance.
- Added pure in-memory EdU and POI normalization APIs.
- Added structured analysis, quantitation, plotting, and safe saving APIs.
- Made failed POI peak detection an error unless the legacy fallback is
  explicitly selected.
- Replaced the active Quarto workflow with a thin package client.
- Isolated optional Python orchestration outside the installed package while
  preserving the Python exporter unchanged.
- Added unit, integration, regression, package-check, and Quarto smoke tests.
- Analyses now retain all supported quantitation for background-subtracted and
  normalized signal representations.
- Added durable analysis and editable figure-bundle RDS APIs.
- Added dynamic, validated presentation overrides and partial condition-color
  mappings that never require templates to know sample labels in advance.
- Added focused complete, pseudocolor, quantitation, cell-cycle, and diagnostic
  Quarto templates with Plotgardener multipanel assembly.
- Changed the automatic condition-color palette from pink to viridis.
