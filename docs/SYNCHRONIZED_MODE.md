# Synchronized mode

`plot_type: synchronized` is a reusable DNA-normalization mode for explicitly
mapped samples. It does not infer gates, alter the target channel, calculate
phase bins, or render an experiment-specific report.

For every configured sample, provide a FlowJo-exported `Single Cells` CSV with
a unique, non-missing `event_index` column. By default, also provide a matching
FlowJo-exported `G1` CSV with the same stable source-event identifiers. The G1
event indices must be an exact subset of the matching Single Cells event
indices; missing, duplicate, or mismatched identifiers stop analysis. This
mode does not attempt to establish membership from measured channel values.

```yaml
plot_type: synchronized
data_dir: path/to/exports
dna_channel: FL2-A
target_channel: FL4-A
target_name: EdU
output_pdf: results/panels.pdf
output_png: results/panels.png
synchronized_dna_strategy: per_sample_g1
synchronized_minimum_g1_events: 100
replicates:
  - label: Replicate 1
    samples:
      - label: Control
        prefix: control_rep1
      - label: Treatment
        prefix: treatment_rep1
```

The defaults expect `<prefix>_single_cells.csv` and `<prefix>_g1.csv`. To use
one asynchronous-reference G1 population for all samples, set
`synchronized_dna_strategy: shared_asynchronous_g1` and provide an explicit
`synchronized_reference_prefix`. This is a scientific-method choice and must
be recorded in the experiment configuration.

The analysis sets `dna_norm` from the approved G1 DNA anchor and copies the
acquired target signal unchanged to `target_raw`, `target_norm`, and
`target_bgsub`. A report, gate definition, sample mapping, and downstream
quantitation remain separate, explicit work.
