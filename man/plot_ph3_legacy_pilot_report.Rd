% Generated manually for the limited-provenance pilot API.
\name{plot_ph3_legacy_pilot_report}
\alias{plot_ph3_legacy_pilot_report}
\title{Plot a limited-provenance legacy pH3 pilot report}
\usage{
plot_ph3_legacy_pilot_report(analysis)
}
\arguments{
\item{analysis}{A completed \code{legacy_csv_pilot_v1} pH3 analysis.}
}
\value{A named list of four outcome ggplot panels, one background-subtracted pH3-versus-DNA
pseudocolor panel per configured pilot sample, and limited-provenance QC tables.}
\description{Produces a clearly labelled pilot-only report from existing legacy
FlowJo CSV exports. It does not claim verified event identity or containment.
Each sample first receives an auditable provisional raw-4N partition used only
to fit a DNA-dependent linear background model. Final cutoff, positivity, C/D
signal values, and displayed pH3 signal use the background-subtracted scale;
the final cutoff is derived only from the matched Untreated control.}
\details{This function is solely for an explicitly configured pilot smoke test.
It must not be used for publication-grade or final analyses.}
