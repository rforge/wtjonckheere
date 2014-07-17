\name{wtJonckheere-package}
\alias{wtJonckheere-package}
\alias{wtJonckheere}
\docType{package}
\title{
  Non-parametric trend test with weighted data points
}
\description{
  The Jonckheere Terpstra test evaluates an trend of data point values between ordered groups.
  This package extends this test procedure with the capability to associate weights
  for individual data points. Weighted data points may be grouped together.
}
\details{
\tabular{ll}{
Package: \tab wtJonckheere\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-07-10\cr
License: \tab GPL-2\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Andreas Recke

Maintainer: Andreas Recke <Andreas.Recke@gmx.de>
}
\references{
    Mack, H. B. (1981)
    K-sample rank tests for umbrella alternatives.
    \emph{J. Amer. Statist. Ass.} \bold{76}, 175--181.

    Jonckheere, A. R. (1954)
    A distribution-free k-sample test against ordered alternatives.
    \emph{Biometrica} \bold{41}, 133--145.

    Terpstra, T. J. (1952)
    The asymptotic normality and consistency of Kendall's test against trend, when ties are present in one ranking.
    \emph{Indagationes Mathematicae} \bold{14}, 327-333.
}
\keyword{trend test}
\keyword{non-parametric}
\seealso{
 \code{\link[clinfun]{jonckheere.test}}
 \code{\link[ENMISC]{wtd.boxplot}}
}
% \examples{
% %% ~~ simple examples of the most important functions ~~
% }
