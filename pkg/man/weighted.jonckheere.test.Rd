\name{weighted.jonckheere.test}
\alias{weighted.jonckheere.test}
\docType{package}
\title{
Weighted Jonckheere Terpstra testing
}
\description{
A permutation based non-parametric trend testing
}
\usage{
weighted.jonckheere.test(x, g, wt, id = factor(g), B=1e4, normalize.wts = FALSE, alternative = c("two.sided", "increasing", "decreasing"))
}
\examples{
\dontrun{
  # define groups
  groups = rep(1:5, each=3)
  values = sapply(groups, function(x) rnorm(1, mean=x, sd=0.5))

  # simulate some weights
  weights = runif(length(groups))

  wttest <- weighted.jonckheere.test(x=values, g=groups, wt = weights, B=1e4, alternative = "two.sided")
}
}


 