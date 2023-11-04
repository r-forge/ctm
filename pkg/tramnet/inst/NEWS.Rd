
\name{NEWS}
\title{NEWS file for the tramnet package}

\section{Changes in Version 0.1-0 (2023-11-04)}{
  \itemize{
    \item{Migrated from \code{mlrMBO} to \code{mlr3mbo} due to discontinued
      support.}
    \item{New formula interface in \code{tramnet.formula()}. Syntax with
      unconditional \code{"tram"} object and design matrix \code{x} can
      still be used; now implemented in \code{tramnet.tram()}.}
    \item{New functions \code{LmNET()}, \code{ColrNET()}, \code{BoxCoxNET()},
      \code{SurvregNET()}, \code{LehmannNET()} and \code{CoxphNET()} using
      the new formula interface.}
    \item{Better handling of out-of-sample evaluation with \code{newdata}.}
  }
}

\section{Changes in Version 0.0-3 (2020-05-14)}{
  \itemize{
    \item{Correct citations in \code{DESCRIPTION}.}
    \item{Correct lexicographical errors in vignette and documentation.}
  }
}

\section{Changes in Version 0.0-2 (2020-03-20)}{
  \itemize{
    \item{\pkg{tramnet} first published on CRAN.}
  }
}
