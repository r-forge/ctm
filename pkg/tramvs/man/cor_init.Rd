% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tramvs.R
\name{cor_init}
\alias{cor_init}
\title{Compute correlation for initializing the active set}
\usage{
cor_init(m0, mb)
}
\arguments{
\item{m0}{\code{modFUN(formula, data)}}

\item{mb}{\code{modFUN(mandatory, data)}}
}
\value{
Vector of correlations for initializing the active set, depends on
    type of model (see e.g. \code{\link[tramvs]{cor_init.default}})
}
\description{
Compute correlation for initializing the active set
}
