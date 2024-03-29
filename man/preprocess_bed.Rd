% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_bed.R
\name{preprocess_bed}
\alias{preprocess_bed}
\title{preprocess_bed preprocess bed file to remove windows of size maxlag for chich
there is no signal or only signal in one sample}
\usage{
preprocess_bed(
  files,
  names = files,
  maxlag = 1e+05,
  prefix = "",
  outdir = tempdir(check = TRUE)
)
}
\arguments{
\item{files}{a vector of bed files path}

\item{names}{(default: files) a list of names for the sample}

\item{maxlag}{(default: 1e5) the size of the windows}

\item{prefix}{(default: "clean_") the prefix of the preprocessed bed file}

\item{outdir}{(default: "tempdir()") path for the preprocessed bed files}
}
\value{
save a list of preprocessed bed files
}
\description{
preprocess_bed preprocess bed file to remove windows of size maxlag for chich
there is no signal or only signal in one sample
}
