
\name{NEWS}
\title{NEWS file for the tram package}

\section{Changes in Version 1.0-0 (2023-07-13)}{
  \itemize{
    \item{Improved, simplified, and generalised implementation of \code{mmlt}.}
    \item{\code{Coxph}, \code{Colr}, and \code{BoxCox} allow unbiased
          estimation of time-to-event distributions under dependent
          censoring.}
    \item{All \code{tram} models can be equipped with frailty terms.}
  }
}

\section{Changes in Version 0.8-2 (2023-03-28)}{
  \itemize{
    \item{Update demo \code{stram}.}
  }
}



\section{Changes in Version 0.8-1 (2023-03-09)}{
  \itemize{
    \item{Add \code{simulate.mmlt}.}
    \item{\code{Survreg} didn't work with constraints, reported by Lucas Kook.}
    \item{\code{mtram} didn't deal properly with interval censoring in all cases.}
    \item{Enable stratification for shift-scale transformation models (still experimental).}
    \item{Fix S3 argument issues.}
    \item{Try to detect negative SEs.}
  }
}

\section{Changes in Version 0.8-0 (2022-11-02)}{
  \itemize{
    \item{Stabilise and generalise \code{mmlt}: Alternative link functions, 
          separate lower triangular matrix storage and computations from
          \code{mmlt}.}
    \item{\code{coef(<mmlt>, type = "Sigma")} and \code{coef(<mmlt>, type = "Corr")}
          were incorrect for more than three dimensions.}
    \item{Improved numerical stability in \code{mtram}.}
  }
}

\section{Changes in Version 0.7-2 (2022-08-05)}{
  \itemize{
    \item{New location-scale functionality, see \code{demo("stram")}.}
    \item{New Laplace link function.}
    \item{New cauchit link function.}
  }
}


\section{Changes in Version 0.7-1 (2022-05-18)}{
  \itemize{
    \item{Minor bug fixed for shift-scale transformation models (still experimental).}
    \item{\code{bread.tram} added for \code{sandwich} to work.}
    \item{Update \code{mtram} vignette; move code to demo.}
  }
}


\section{Changes in Version 0.7-0 (2022-04-01)}{
  \itemize{
    \item{Add shift-scale transformation models (experimental).}
    \item{Add methods for probabilistic indices and overlap coefficients + friends.}
    \item{Export \code{Gradient} and \code{Hessian} methods.}
    \item{Simplify \code{score_test} and \code{profile}: Use \code{fixed}
          instead of \code{offset}.}
    \item{Models with standardisation now available in \code{mtram} for exact
          continuous outcomes and link functions other than probit.}
  }
}

\section{Changes in Version 0.6-4 (2022-01-14)}{
  \itemize{
    \item{Remove arguments ignored by \code{mlt}.}
    \item{Move mtram vignette code to demo as it took too long 
          on Windows for CRAN.}
  }
}

\section{Changes in Version 0.6-3 (2022-01-04)}{
  \itemize{
    \item{Deal with stratum - covariate interactions,
          overparameterised coefficients are constraint to zero.}
  }
}

\section{Changes in Version 0.6-2 (2021-10-19)}{
  \itemize{
    \item{Make sure to use the Schur complement when computing
          variances in for linear transformation models in 
          \code{score_test} and \code{perm_test}.}
    \item{Fix LaTeX problems.}
  }
}

\section{Changes in Version 0.6-1 (2021-09-16)}{
  \itemize{
    \item{Add \code{dofit} argument to \code{mmlt}.}
    \item{Robustify \code{vcov} for \code{mmlt} objects.}
    \item{Add experimental \code{score_test} and \code{perm_test} methods
          for \code{glm}.}
  }
}


\section{Changes in Version 0.6-0 (2021-03-08)}{
  \itemize{
    \item{Add \code{simulate} method.}
    \item{Enable nonparametric likelihood estimation by \code{Polr}.}
    \item{New argument \code{sparse_nlevels} to \code{tram}, requesting
          a model \code{Matrix} for ordered factors with at least
          \code{sparse_nlevels} levels.}
  }
}


\section{Changes in Version 0.5-2 (2020-10-23)}{
  \itemize{
    \item{Add experimental \code{score_test} and \code{perm_test} for
          \code{survival::coxph}.}
    \item{Update simulation code.}
  }
}


\section{Changes in Version 0.5-1 (2020-09-03)}{
  \itemize{
    \item{Robustify computation of score intervals.}
    \item{Add Spearman correlation to \code{mmlt} example.}
  }
}

\section{Changes in Version 0.5-0 (2020-08-27)}{
  \itemize{
    \item{\code{Lehmann()} parameter is now labelled as proportional reverse time hazard ratio.}
    \item{Add \code{mmlt} for multivariate conditional transformation models
          (continuous responses only for now).}
    \item{Reproducibility material for <arxiv:1906.03151> added to package.}
  }
}


\section{Changes in Version 0.4-0 (2020-05-11)}{
  \itemize{
    \item{Force support to start at time zero and constrain cumulative
          hazard to zero for time zero in \code{Aareg}.}
  }
}


\section{Changes in Version 0.4-0 (2020-05-11)}{
  \itemize{
    \item{Add \code{Aareg} for Aalen additive hazard models.}
  }
}


\section{Changes in Version 0.3-3 (2020-01-16)}{
  \itemize{
    \item{\code{score_test} generic and method added for computing 
          score tests and corresponding confidence intervals.}
    \item{\code{perm_test} generic and method added for computing 
          permutation score tests.}
    \item{\code{PI} converts log-odds ratios to probabilistic indices and 
          vice versa.}
  }
}


\section{Changes in Version 0.3-2 (2020-01-08)}{
  \itemize{
    \item{\code{profile} method added for computing 
          profile likelihood confidence intervals.}
    \item{Allow additional linear inequality constraints for regression
          coefficients.}
  }
}

\section{Changes in Version 0.3-1 (2019-11-05)}{
  \subsection{Bugfixes}{
    \itemize{
      \item{Fix problem in C header file.}
      \item{Fix problems with responses of class \code{R}, spotted by Balint Tamasi.}
    }
  }
}

\section{Changes in Version 0.3-2 (2019-12-02)}{
  \itemize{
    \item{Add new demo.}
  }
  \subsection{Bugfixes}{
    \itemize{
      \item{Deal with \code{model_only} in \code{Survreg}, spotted by Balint Tamasi.}
    }
  }
}


\section{Changes in Version 0.3-0 (2019-10-09)}{
  \subsection{Bugfixes}{
    \itemize{
      \item{\code{y | s1:s2 ~ x} lead to rank-deficit models, the intercept is
            now removed.}
      \item{Not really a bug, but time-consuming examples were removed
            (or replaced by static output) in package vignettes to make
            CRAN happy.}
      \item{\code{as.mlt} returned an object of class \code{mlt} but not
            \code{mlt_fit}.}
    }
  }
  \subsection{New Features}{
    \itemize{
      \item{New marginally interpretable models for clustered data
            in \code{mtram}, see \code{mtram} vignette. This is still
            experimental and very rough, so please be careful.}
      \item{Compute \code{vcov.tram} via Schur complement when possible.}
  }
  }
}

\section{Changes in Version 0.2-6 (2019-06-20)}{
  \itemize{
    \item{Add regression test for count data and \code{log_first = TRUE}.}
  }
}


\section{Changes in Version 0.2-5 (2019-02-25)}{
  \itemize{
    \item{\code{tram} objects are now in line with \code{multcomp::glht},
          allowing post-hoc multiple comparisons.}
  }
}

\section{Changes in Version 0.2-4 (2019-01-22)}{
  \itemize{
    \item{\code{log_first = TRUE} is now accepted by \code{tram} and all 
          model functions calling \code{tram}. This is still experimental.}
    \item{Silence regression tests on Windows.}
  }
}


\section{Changes in Version 0.2-3 (2018-09-17)}{
  \itemize{
    \item{New \code{extrapolate = TRUE} option for all models constraints
          the second derivate of the Bernstein transformation function to be
          zero at \code{support}. This helps to get rid of possible
          discontinuous densities at \code{support}.}
    \item{Add \code{Lehmann} for fitting shift transformation models those
          linear predictor can be interpreted as log-Lehmann-alternative.}
    \item{Add \code{method = "loglog"} to \code{Polr}.}
    \item{New \code{update} method.}
  }
}

\section{Changes in Version 0.2-2 (2018-06-22)}{
  \itemize{
    \item{Use default \code{bounds = c(0, Inf)} for \code{Surv} objects (for example
          in \code{Coxph}.}
    \item{Update to \pkg{mlt} 1.0-2.}
  }
}


\section{Changes in Version 0.2-1 (2018-04-20)}{
  \itemize{
    \item{Likelihood ratio tests with fixed parameters gave an error.}
    \item{Remove strata from linear predictor (experimental).}
    \item{Try to implement a better guess for the \code{add} argument
          later defining the axes corresponding to the response
          in plots.}
    \item{All models accept \code{add} and \code{bounds} arguments.}
  }
}

\section{Changes in Version 0.2-0 (2018-02-25)}{
  \itemize{
    \item{\pkg{tram} first published on CRAN.}
  }
}
