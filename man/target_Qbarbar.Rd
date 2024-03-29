% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/target.R
\name{target_Qbarbar}
\alias{target_Qbarbar}
\title{Target the iterated regressions}
\usage{
target_Qbarbar(
  Qbar,
  Qbarbar,
  Y,
  A,
  a,
  a_star,
  gn,
  M1,
  M2,
  all_mediator_values,
  max_iter,
  tol = 1/(sqrt(length(Y)) * log(length(Y)))
)
}
\arguments{
\item{Qbar}{Outcome regression estimates}

\item{Qbarbar}{Iterated mean estimates}

\item{Y}{A vector of continuous or binary outcomes.}

\item{A}{A vector of binary treatment assignment (assumed to be equal to 0 or
1).}

\item{a}{The label for the treatment. The effects estimates returned pertain
to estimation of interventional effects of \code{a} versus \code{a_star}.}

\item{a_star}{The label for the treatment. The effects estimates returned pertain
to estimation of interventional effects of \code{a} versus \code{a_star}.}

\item{gn}{Power users may wish to pass in their own properly formatted list of the
propensity score so that
nuisance parameters can be fitted outside of \code{intermed}.}

\item{M1}{A \code{vector} of mediators.}

\item{M2}{A \code{vector} of mediators.}

\item{all_mediator_values}{All combinations of M1 and M2}

\item{max_iter}{The maximum number of iterations for the TMLE}

\item{tol}{The tolerance for stopping the iterative targeting procedure.}
}
\value{
A list containing all outcome regression evaluations needed
for downstream calculations.
}
\description{
Target the iterated regressions
}
