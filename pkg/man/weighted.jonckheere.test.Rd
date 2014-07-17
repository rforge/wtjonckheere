\name{weighted.jonckheere.test}
\alias{weighted.jonckheere.test}
\alias{print.weighted_jt_statistic}
\alias{plot.weighted_jt_statistic}
\alias{summary.weighted_jt_statistic}
\docType{package}
\title{
Weighted Jonckheere Terpstra testing
}
\description{
A permutation based non-parametric trend testing
}
\usage{
weighted.jonckheere.test(x, g, wt, id = factor(x), B=1e4, normalize.wts = FALSE, alternative = c("two.sided", "increasing", "decreasing"), grouping = c("none","by.groups","by.values"))
}
\arguments{
\item{x}{ a numeric vector containing values }
\item{g}{ an ordered factor or integer vector containing the groups }
\item{wt}{ a numeric vector the probability weight for data point }
\item{id}{ a factor that associates individuals to data points (see details) }
\item{B}{ the number of permutations }
\item{normalize.wts}{ whether the weights should be normalized to 1.0 for data points associated with single \code{id}s}
\item{alternative}{ whether the trend is increasing, decreasing or changing into any direction }
\item{grouping}{ determines the way how permutation takes \code{id} into account (see details) }
}
\details{

  The null hypothesis of a equality between the groups \code{g} concerning the values \code{x} is tested.

  The permutation test implented here is based upon the regular Jonckheere Terpstra test, but incorporates
  weights \code{wt} and the possiblity group data points by individuals \code{id}. When data points are grouped by
  individuals, either the x values (\code{grouping == "by.groups"}) or the groups (\code{by.values}) are randomized
  by the permutation process. if \code{grouping == "none"} is selected, the permutation test will randomize both.
  In all cases, the weights (\code{wt}) remain unrandomized.
  
  The weighted Jonckheere Terpstra statistic \eqn{JT} is calculated for \eqn{n} data points
  \deqn{JT = \sum_{i=1}^n \sum{j=1}^n w_i w_j U_{ij}}
  with \eqn{i \leq n} and \eqn{j \leq n} and \eqn{i, j, n \in \mathcal{N}}, weights \eqn{w}, values \eqn{x} and groups \eqn{g}.
  \deqn{U_{ij} = \begin{cases} 1 \mbox{if } g_i < g_j \mbox{and} x_i < x_j \\
                               0.5 \mbox{if } g_i < g_j \mbox{and} x_i = x_j \end{cases}}

  In case of all weights identical to 1, this statistic is equal to the classical JT statistic.

  The plot method for this function provides an easy overview of the permutation samples.
}
\value{
An object of class \code{weighted_jt_statistic}, basically a list including elements
\item{statistic}{the reference Jonckheere Terpstra statistic}
\item{permutations}{a vector of Jonckheere Terpstra statistics of permutations}
\item{alternative}{the test setting - as determined by input parameters}
\item{mean.permutation}{the mean of Jonckheere Terpstra statistics of permutations}
\item{std.permutation}{the standard deviation of Jonckheere Terpstra statistics of permutations}
\item{n.permutations}{the number of permutations - as determined by input parameters}
\item{est.p.value}{the p value determined from a normal distribution approximated from permutations}
\item{exact.p.value}{the p value determined from the permutated samples}
}

\author{Andreas Recke}
\examples{
  # load a single data set
  data(mrna)
  wttest <- with(mrna, weighted.jonckheere.test(x=mRNA.levels, g=gene.copies, id=id, wt = weights, B=1e4, alternative = "two.sided", grouping = "by.groups"))
  summary(wttest)
  # and plot
  plot(wttest)
}



 