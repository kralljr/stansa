% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-res.R
\name{plotstan}
\alias{plotstan}
\title{plotstan}
\usage{
plotstan(
  typesim,
  stanres,
  dirname = NULL,
  filename = NULL,
  by = 5,
  prof = prof,
  meansd = meansd,
  pdf = F,
  hten = 150,
  wden = 20,
  htrh = 100,
  htbi = 10,
  wdbi = 10
)
}
\arguments{
\item{typesim}{Which simulation to run.  Currently, local1 - local4, ambient.}

\item{stanres}{Output from runstan}

\item{dirname}{Folder name where to save plots}

\item{filename}{Base filename or typesim if NULL}

\item{by}{Skips for traceplot, biasplot}

\item{prof}{Profile data, defaults to prof.}

\item{meansd}{Matrix with columns corresponding to source, type, mean, sd, defaults to meansd.}

\item{pdf}{Whether to save pdfs of energy, rhat, bias}

\item{hten}{Height in inches of energy pdf}

\item{wden}{Width in inches of bias pdf}

\item{htrh}{Height in inches of rhat pdf}

\item{htbi}{Height in inches of bias pdf}
}
\description{
\code{plotstan} Create plots for assessing stan diagnostic
}
