# facspseudocolor 0.0.0.9000

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
