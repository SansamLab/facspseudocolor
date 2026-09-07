# Figure 6 EdU NMPP1 CSV-only pilot setup

**Status:** Real-data model validation complete; HTML visual review pending a
Quarto-enabled environment.

## Scope

The approved candidate is the documented Figure 6 EdU NMPP1 experiment:

`/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/Figure6_FlowJoAndFcsFiles/facspseudocolor_analysis_EdUNmpp1Treatment/config.yml`

The original configuration declares two independent biological-replicate runs
(2025-07-28 and 2025-07-30), each containing Untreated, 5 h auxin, 3 h
1NM-PP1, and 5 h auxin + 3 h 1NM-PP1 samples. The pilot preserves those eight
explicit labels, prefixes, FCS source names, reference labels, channel names,
DNA normalization anchor, EdU background-subtracted quantitative settings, and
display settings. No mapping was inferred from file order.

The new pilot configuration is outside the historical experiment:

`/Users/sansamc/OMRF Dropbox/Chris Sansam/Synchd/Organized Dropbox/Projects and Data/Projects/Cowork_Projects/2026_Analyze_FACS_Data/test_reports/edu_figure6_nmpp1_csv_pilot_2026-09-07/config.yml`

It refers directly to the historical `csv` directory through an explicit
absolute path, specifies `flowjo.rebuild: false`, and sends any configured PDF
or PNG output paths to the pilot directory. The Quarto render command builds
and writes only below a newly created external temporary directory.

## What this pilot does and does not establish

It exercises the newly merged every-sample EdU DNA-versus-EdU pseudocolor
report against a compact, documented real experiment. It does not alter the
historical CSV exports or FlowJo/FCS inputs, regenerate exports, revise any
mapping, or replace the established quantitative report workflow. Its purpose
is visual inspection of the EdU panels plus display-offset, suppression, and
provenance tables.

## Inputs and outputs

Experimental source contents inspected while preparing this setup: the
original YAML configuration and file names only. The subsequent approved pilot
read the 24 explicitly named existing CSV exports to construct its analysis and
report model. No FCS, WSP, RDS, PDF, PNG, or historical HTML content was read;
no experimental input was modified.

Files generated:

- `test_reports/edu_figure6_nmpp1_csv_pilot_2026-09-07/config.yml`
- `test_reports/edu_figure6_nmpp1_csv_pilot_2026-09-07/README.md`
- this record.

No generated repository artifact was created. The approved pilot created its
package build and installation only beneath an external temporary directory;
it constructed all eight named report panels with status `available`. The
sandbox did not have Quarto, so it did not create HTML there. No historical
experiment artifact was written or modified.

## Scientific scope

Sample mappings, exclusions, gates, thresholds, normalization, statistics, and
biological claims changed: **none**. The user approved only a safe CSV-only
pilot configuration using the existing documented mapping.

## Verification

The setup task itself did not run code. Under the temporary execution
exception, the final configuration correction was then verified by focused
configuration tests (PASS 30), the full package suite (PASS 1,753), and
`R CMD check --no-manual` (Status OK). The approved real-data model
construction completed with all eight panels available. Its sole warning was
the known legacy EdU table-alias deprecation; no new warning or fallback
occurred. The exact Quarto render command and stop conditions remain in the
pilot `README.md`.

## Assumptions and unresolved points

This setup assumes that all 24 declared CSV exports remain present and have the
documented channel names and membership expected by the package. Those contents
were intentionally not inspected. A missing input or validation failure must
stop the pilot; it must not select another file or silently omit a sample.
