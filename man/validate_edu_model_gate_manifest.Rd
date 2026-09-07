\name{validate_edu_model_gate_manifest}
\alias{validate_edu_model_gate_manifest}
\title{Validate provenance for experimental model-derived EdU gate inputs}
\usage{validate_edu_model_gate_manifest(analysis, manifest_path)}
\arguments{
\item{analysis}{A completed \code{facs_analysis} in EdU mode.}
\item{manifest_path}{Explicit path to \code{model-gating-manifest.json}.}
}
\value{An invisible data frame describing the verified input artifacts.}
\description{Binds an EdU analysis to the non-production model-gating manifest
and verifies the path and SHA-256 digest of every population input.}
