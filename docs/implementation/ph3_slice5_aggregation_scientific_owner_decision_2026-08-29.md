# pH3 Slice 5 aggregation — scientific-owner decision record

**Decision status:** CONFIRMED

**Confirmation date:** 2026-08-29

**Scientific owner:** Chris Sansam

## Scope

This record resolves the mandatory contract gate encountered before pH3 Slice
5 implementation. It supplements the existing pH3 scientific-owner approval,
implementation plan/test specification, and merged Slice 1–4 implementation
records. It does not reopen or revise the event-level, direct-identity,
containment, positivity, normalization, phase, interval, acquisition-metric, or
boundary-sensitivity methods.

## Authoritative source values

Slice 5 consumes only the four retained Slice 4 acquisition tables. Numerical
aggregation is performed only for:

- `ph3_metrics_acquisition`;
- `ph3_phase_prevalence`.

The following remain acquisition-level QC and are not numerically aggregated:

- `ph3_event_eligibility_qc`;
- `ph3_4n_boundary_sensitivity_qc`.

All four Slice 4 acquisition tables remain attached unchanged. Slice 5 never
recomputes event classification, positivity, eligibility, configured phase,
sub-4N/4N ownership, acquisition percentages, or sensitivity predicates.

## Experimental unit and aggregation rules

### Technical acquisitions within a biological replicate

Groups are defined only by the explicit `condition` × `replicate` manifest
mapping. Technical-acquisition relationships must never be inferred from row
order, filenames, prefixes, or adjacency.

For each metric or configured-phase row, the biological-replicate value is the
unweighted arithmetic mean of finite Slice 4 technical-acquisition percentage
values within the explicit condition × replicate group.

The required method identifier is:

`unweighted_technical_percentage_mean_within_biorep_condition_v1`

Numerators and denominators are not pooled or summed to create this value.
Technical acquisitions do not increase the biological-replicate sample size.

### Biological replicates within a condition

For each metric or configured-phase row, the condition mean is the unweighted
arithmetic mean of finite biological-replicate means. Every biological
replicate has equal weight regardless of its number of technical acquisitions.

The required method identifier is:

`unweighted_biological_replicate_mean_within_condition_v1`

Condition summaries are descriptive only. They report the mean, biological-
replicate SD, biological-replicate SEM, and biological-replicate counts. Slice
5 adds no hypothesis test, p-value, confidence interval, multiple-testing
procedure, or biological claim.

## Approved output tables and attachment

Append exactly these four tables beneath `analysis$quantitation`, after the
four retained Slice 4 tables:

1. `ph3_metrics_biological_replicate`;
2. `ph3_phase_prevalence_biological_replicate`;
3. `ph3_metrics_condition_summary`;
4. `ph3_phase_prevalence_condition_summary`.

No compatibility aliases or nested alternative copies are approved.

All four Slice 5 tables use:

`output_schema_version = ph3-1.0.0`

## Exact biological-replicate table schemas

### `ph3_metrics_biological_replicate`

Columns, in order:

1. `analysis_id`
2. `condition`
3. `condition_index`
4. `replicate`
5. `replicate_index`
6. `aggregation_level`
7. `metric_id`
8. `population_id`
9. `interval_lower`
10. `interval_upper`
11. `lower_inclusive`
12. `upper_inclusive`
13. `value_percent`
14. `technical_acquisition_count`
15. `finite_technical_acquisition_count`
16. `undefined_technical_acquisition_count`
17. `result_status`
18. `source_acquisition_ids`
19. `source_export_operation_ids`
20. `source_input_manifest_keys`
21. `classification_schema_version`
22. `output_schema_version`
23. `aggregation_method_id`
24. `positivity_method_id`
25. `eligibility_method_id`
26. `interval_method_id`
27. `four_n_method_id`
28. `sub_four_n_method_id`
29. `containment_method_id`
30. `config_digest`

### `ph3_phase_prevalence_biological_replicate`

Columns, in order:

1. `analysis_id`
2. `condition`
3. `condition_index`
4. `replicate`
5. `replicate_index`
6. `aggregation_level`
7. `metric_id`
8. `phase_id`
9. `phase_index`
10. `interval_lower`
11. `interval_upper`
12. `lower_inclusive`
13. `upper_inclusive`
14. `value_percent`
15. `technical_acquisition_count`
16. `finite_technical_acquisition_count`
17. `undefined_technical_acquisition_count`
18. `result_status`
19. `source_acquisition_ids`
20. `source_export_operation_ids`
21. `source_input_manifest_keys`
22. `classification_schema_version`
23. `output_schema_version`
24. `aggregation_method_id`
25. `positivity_method_id`
26. `eligibility_method_id`
27. `interval_method_id`
28. `four_n_method_id`
29. `sub_four_n_method_id`
30. `containment_method_id`
31. `config_digest`

For both biological-replicate tables:

- `aggregation_level` is `biological_replicate`;
- `aggregation_method_id` is
  `unweighted_technical_percentage_mean_within_biorep_condition_v1`;
- `source_acquisition_ids`, `source_export_operation_ids`, and
  `source_input_manifest_keys` are deterministic canonical JSON arrays of
  unique strings sorted lexicographically;
- the retained acquisition tables are the lossless record of each source
  value and status; Slice 5 does not duplicate source values into a JSON field;
- acquisition numerator/denominator counts are not copied, summed, pooled, or
  represented as biological-replicate counts;
- no technical-acquisition SD or SEM column is added. Technical heterogeneity
  remains auditable from the retained acquisition rows.

## Exact condition-summary table schemas

### `ph3_metrics_condition_summary`

Columns, in order:

1. `analysis_id`
2. `condition`
3. `condition_index`
4. `aggregation_level`
5. `metric_id`
6. `population_id`
7. `interval_lower`
8. `interval_upper`
9. `lower_inclusive`
10. `upper_inclusive`
11. `mean_percent`
12. `sd_percent`
13. `sem_percent`
14. `biological_replicate_count`
15. `finite_biological_replicate_count`
16. `undefined_biological_replicate_count`
17. `result_status`
18. `source_replicate_ids`
19. `classification_schema_version`
20. `output_schema_version`
21. `aggregation_method_id`
22. `positivity_method_id`
23. `eligibility_method_id`
24. `interval_method_id`
25. `four_n_method_id`
26. `sub_four_n_method_id`
27. `containment_method_id`
28. `config_digest`

### `ph3_phase_prevalence_condition_summary`

Columns, in order:

1. `analysis_id`
2. `condition`
3. `condition_index`
4. `aggregation_level`
5. `metric_id`
6. `phase_id`
7. `phase_index`
8. `interval_lower`
9. `interval_upper`
10. `lower_inclusive`
11. `upper_inclusive`
12. `mean_percent`
13. `sd_percent`
14. `sem_percent`
15. `biological_replicate_count`
16. `finite_biological_replicate_count`
17. `undefined_biological_replicate_count`
18. `result_status`
19. `source_replicate_ids`
20. `classification_schema_version`
21. `output_schema_version`
22. `aggregation_method_id`
23. `positivity_method_id`
24. `eligibility_method_id`
25. `interval_method_id`
26. `four_n_method_id`
27. `sub_four_n_method_id`
28. `containment_method_id`
29. `config_digest`

For both condition-summary tables:

- `aggregation_level` is `condition`;
- `aggregation_method_id` is
  `unweighted_biological_replicate_mean_within_condition_v1`;
- `source_replicate_ids` is a deterministic canonical JSON array of unique
  replicate strings sorted lexicographically;
- biological-replicate numerator/denominator counts are not created;
- `sd_percent` is the sample SD across finite biological-replicate means;
- `sem_percent = sd_percent / sqrt(finite_biological_replicate_count)` when at
  least two finite biological replicates exist;
- SD and SEM are `NA_real_` when fewer than two finite biological replicates
  exist.

## Undefined and partial-result behavior

At either aggregation level:

- all source values finite: calculate the approved unweighted mean and use
  `result_status = ok`;
- at least one finite and at least one undefined source value: calculate the
  mean from finite values only and use
  `result_status = ok_partial_undefined`;
- no finite source values: return `NA_real_` and use
  `result_status = undefined_no_finite_values`.

At condition level, exactly one finite biological replicate is a valid
descriptive mean. Its SD and SEM are `NA_real_`. When it is the only declared
biological replicate and no replicate is undefined, use:

`result_status = ok_single_biological_replicate`

If one finite biological replicate coexists with one or more undefined
biological replicates, the partial-undefined status takes precedence:

`result_status = ok_partial_undefined`

Numeric zero remains zero. It is not treated as missing or undefined.

## Counts, completeness, exclusions, and minimum replication

- Acquisition numerator and denominator fields remain only in Slice 4
  acquisition tables.
- Aggregate tables report only total, finite, and undefined source-unit counts.
- No new exclusion mechanism is added.
- Every acquisition explicitly included by the active sample manifest must
  contribute exactly one required row for every approved metric or phase.
- Missing, duplicated, extra, mixed-analysis, mixed-method, mixed-schema, or
  provenance-inconsistent rows are fatal.
- One biological replicate is sufficient for a descriptive condition mean but
  not for SD, SEM, inferential statistics, or a replication claim.
- Zero declared acquisitions or zero declared biological replicates for an
  expected group is fatal rather than an empty or fabricated result.

## Deterministic ordering

- `ph3_metrics_biological_replicate`: `condition_index`, `replicate_index`, then
  the fixed Slice 4 five-metric order.
- `ph3_phase_prevalence_biological_replicate`: `condition_index`,
  `replicate_index`, then `phase_index`.
- `ph3_metrics_condition_summary`: `condition_index`, then the fixed Slice 4
  five-metric order.
- `ph3_phase_prevalence_condition_summary`: `condition_index`, then
  `phase_index`.
- Canonical source-ID arrays use lexicographic ordering only for provenance
  serialization. They never establish sample, condition, replicate, or
  technical-acquisition membership.

## Provenance and fail-closed requirements

Every aggregation group must have exactly consistent values for all retained
method, schema, configuration, and analysis provenance fields. Source
acquisition/operation/manifest identifiers must reconcile to the active Slice
4 tables and explicit manifest. Mixed or missing provenance is fatal. No
filename, row-order, nearest-file, sequential-identity, partial-result, or
legacy fallback is permitted.

## Explicitly deferred

The following remain outside Slice 5:

- inferential statistics;
- normalization to a reference condition;
- CSV/RDS export;
- plots and report rendering;
- geometry;
- GUI/configurator behavior;
- compatibility aliases;
- alteration of legacy pH3, EdU, or POI behavior;
- any Slice 6 or later output.

## Confirmation

The scientific owner explicitly confirmed the proposed contract with:

> CONFIRM

This record is therefore authoritative for Slice 5 implementation. Any change
to these scientific semantics requires a new explicit scientific-owner
decision; passing tests alone does not authorize a method change.
