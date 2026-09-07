% Generated manually for the EdU Slice 1 output-contract API.
\name{build_edu_pseudocolor_output_contract}
\alias{build_edu_pseudocolor_output_contract}
\title{Build every-sample EdU pseudocolor panels with display-only offsets}
\usage{
build_edu_pseudocolor_output_contract(analysis)
}
\arguments{
\item{analysis}{A completed EdU \code{facs_analysis} object.}
}
\value{A list of individually identified \code{ggplot} panels and per-sample
display-offset and panel-QC tables.}
\description{Uses the exact EdU-negative event membership retained from each
completed background fit to calculate \code{median(raw) - median(corrected)}
for that panel only. Analytical values, classifications, and exports are not
modified. A panel is explicitly suppressed if its identity, membership proof,
offset, or display coordinate cannot be validated.}
