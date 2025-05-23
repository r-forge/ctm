
\name{NEWS}
\title{NEWS file for the basefun package}

\section{Changes in Version 1.2-3 (2025-05-09)}{
  \itemize{
    \item{\code{rep} to \code{rep_len}.}
  }
}

\section{Changes in Version 1.2-2 (2025-01-29)}{
  \itemize{
    \item{Switch from \pkg{coneproj} to \pkg{quadprog} in \code{?Bernstein_basis}.}
  }
}

\section{Changes in Version 1.2-1 (2025-01-17)}{
  \itemize{
    \item{Allow to keep unused levels in \code{as.basis.formula}.}
  }
}


\section{Changes in Version 1.2-0 (2024-05-21)}{
  \itemize{
    \item{New \code{cyclic_basis()} for modelling seasonal effects.}
  }
}

\section{Changes in Version 1.1-4 (2023-05-16)}{
  \itemize{
    \item{\code{Bernstein_basis(..., ui = c("increasing", "zerointegral")} came with incorrect constraint matrix.}
  }
}

\section{Changes in Version 1.1-3 (2022-05-25)}{
  \itemize{
    \item{\code{Bernstein_basis(..., ui = "zerointegral")} came with incorrect constraint matrix.}
  }
}

\section{Changes in Version 1.1-2 (2022-01-14)}{
  \itemize{
    \item{Allow coefficients in matrix form.}
  }
}


\section{Changes in Version 1.1-1 (2021-11-12)}{
  \itemize{
    \item{Add convexity and concavity constraints to 
          polynomials in Bernstein form.}
    \item{Add toy example demonstrating constrainted 
          curve estimation to \code{?Bernstein_basis}.}
  }
}

\section{Changes in Version 1.1-0 (2021-03-09)}{
  \itemize{
    \item{More infrastructure for \code{as.basis} to return 
          a \code{Matrix} object for ordered factors.}
  }
}

\section{Changes in Version 1.0-7 (2020-02-13)}{
  \itemize{
    \item{Update citation info.}
  }
}

\section{Changes in Version 1.0-6 (2019-10-05)}{
  \itemize{
    \item{Change name of parameter in \code{intercept\_basis} to 
          \code{(_Intercept_)}; this helps to detect 
          intercepts originating from other calls to \code{model.matrix}.}.
  }
}


\section{Changes in Version 1.0-5 (2019-06-20)}{
  \itemize{
    \item{Make sure \code{bounds[1]} is positive if \code{log_first = TRUE}.}
    \item{Correct behaviour outside support when \code{log_first = TRUE}.}
  }
}


\section{Changes in Version 1.0-4 (2019-02-25)}{
  \itemize{
    \item{\code{model.matrix} has \code{constrasts} attribute where appropriate.}
  }
}

\section{Changes in Version 1.0-3 (2019-01-16)}{
  \itemize{
    \item{Bernstein(log()) is possible with new \code{log_first} argument.}
    \item{Fix > length 1 issues.}
  }
}

\section{Changes in Version 1.0-2 (2018-09-17)}{
  \itemize{
    \item{Constrain Bernstein to zero second derivative at the support when
          linear extrapolation shall be used. This avoids ugly bumps in the
          corresponding densities in \pkg{mlt}.}
  }
}


\section{Changes in Version 1.0-1 (2018-03-27)}{
  \itemize{
    \item{Replace \code{rBind} by \code{rbind} as requested by MM.}
  }
}


\section{Changes in Version 1.0-0 (2018-03-07)}{
  \itemize{
    \item{A paper describing version 1.0-0 of the \pkg{mlt}, \pkg{basefun}, and
          \pkg{variables} packages was accepted for publication in the Journal of
          Statistical Software 2018-03-05.}
    \item{Documentation updates.}
  }
}


\section{Changes in Version 0.0-40 (2018-02-23)}{
  \itemize{
    \item{\code{nparm} sometimes generated data with sanitised variable
          names, we don't want this to happen.}
    \item{Allow positive increasing Bernstein polynomials.}
  }
}

\section{Changes in Version 0.0-39 (2017-12-18)}{
  \itemize{
    \item{Add \code{intercept_basis}, ie an intercept only model matrix.}
  }
}


\section{Changes in Version 0.0-38 (2017-05-15)}{
  \itemize{
    \item{\code{model.matrix.cbind_bases} can now deal with 
          lists as \code{newdata} argument (given appropriate 
          dimensions in \code{dim}).}
  }
}



\section{Changes in Version 0.0-37 (2017-04-04)}{
  \itemize{
    \item{Deal with integer variables when \code{sumconstr = TRUE}.}
    \item{Remove unnecessary warning when \code{sumconstr = TRUE} and 
          \code{negative = TRUE}.}
  }
}


\section{Changes in Version 0.0-36 (2017-03-01)}{
  \itemize{
    \item{Register C routines.}
  }
}


\section{Changes in Version 0.0-35 (2017-02-14)}{
  \itemize{
    \item{Fix bug in \code{predict} when \code{newdata} was a data frame.}
  }
}


\section{Changes in Version 0.0-34 (2016-12-14)}{
  \itemize{
    \item{Fix bug in \code{predict} when \code{newdata} was a list. The
          ordering of the resulting array may have been incorrect under
          rare circumstances.}
  }
}

\section{Changes in Version 0.0-33 (2016-09-12)}{
  \itemize{
    \item{Approximate Bernstein polynomial f outside its \code{support} (l, u)
          by a linear function f(u + x) = f(u) + f'(u) * (u - x), now also
          for \code{deriv > 0} (\code{maxderiv} argument in 
          \code{model.matrix.Bernstein_basis}. The default (linear extrapolation
          for \code{deriv = 0} only) did not change.}
  }
}

\section{Changes in Version 0.0-32 (2016-05-31)}{
  \itemize{
    \item{Contraints in \code{log_basis} were potentially wrong.}
  }
}


\section{Changes in Version 0.0-31 (2016-04-25)}{
  \itemize{
    \item{\code{log_basis} has better checks for positivity of
          the variable to be log-transformed.}
  }
}

\section{Changes in Version 0.0-30 (2016-02-19)}{
  \itemize{
    \item{The \pkg{basefun} package was first 
          published on CRAN.}
  }
}
