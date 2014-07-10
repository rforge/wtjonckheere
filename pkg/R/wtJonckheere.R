# The statistic is calculated under the precondition that for each individual id the value x is known but the group g is not clear
# therefore. This is important for 
# It is not checked whether 

weighted.jonckheere.test <- function(x, g, wt, id = factor(g), B=1e4, normalize.wts = FALSE, alternative = c("two.sided", "increasing", "decreasing"))
{
  # check argument
  alternative = match.arg(alternative, c("two.sided", "increasing", "decreasing"))

  # Clean factor "id"
  id <- factor(as.character(id))
  
  # do a check: x must be identical for all groups
  check.v = tapply(x, id, max)
  if (any( x != check.v[id])) stop("Error: x values must be identical for each id individual!\n");
  
  if (normalize.wts) {
    wt = wt / tapply(wt, id, sum)[id]
  }
    
  # First, we calculate the statistic under the given condition
  basis = .Call("rcpp_wt_jonckheere", values = as.double(x), group = as.integer(g), weights = as.double(wt), PACKAGE="wtJonckheere")
  
  # Second, we calculate multiple permutations. For this purpose, the values for x are exchanged between the ids
  perm_x_values = tapply(x, id, mean)
  
  permutationen <- sapply(1:B, function(dum) {
    perm_x_new_values = perm_x_values[order(runif(length(perm_x_values)))]
    names(perm_x_new_values) <- names(perm_x_values)
    return(.Call("rcpp_wt_jonckheere", values = as.double(perm_x_new_values[id]), group = as.integer(g), weights = as.double(wt), PACKAGE="wtJonckheere"))
  })
  
  mittelwert = mean(permutationen)
  std = sd(permutationen)
  pv = pnorm(basis, mittelwert, std)
  
  if (alternative == "increasing")  {
    exakt.pv = sum(basis > permutationen)/B
    pv = 1-pv
  }
  if (alternative == "decreasing") {
    exakt.pv = sum(basis < permutationen)/B    
  }
  if (alternative == "two.sided")
  {
    if (pv > .5) pv = 2*(1-pv) else pv = 2*pv;
    exakt.pv = sum(basis > permutationen)/B  
    if (exakt.pv > .5) exakt.pv = 2*(1-exakt.pv) else exakt.pv = 2*exakt.pv
  }
  
  result = list(statistic = basis, permutationen=permutationen, alternative = alternative, 
                mean.permutation = mittelwert, std.permutation = std, n.permutation=B,
                est.p.value = pv, exact.p.value = exakt.pv);
  class(result) = "weighted_jt_statistic"
  return(result)
}

# And a printing function
summary.weighted_jt_statistic <- function(x, ...)
{
  cat("\nWeighted Jonckheere Terpstra Testing\n")
  cat("JT statistic of data = ", x$statistic,"\n")
  cat("Summary of ", x$n.permutation, " permutations. Average JT statistic =",x$mean.permutation, " +/- std. dev. ", x$std.permutation, "\n")
  qq = quantiles(x$permutationen, probs=c(0.025,0.5,0.975))
  cat("2.5% quantile = ", qq[1], " median= ", qq[2], " 97.5% quantile = " ,qq[3],"\n")
  cat("Alternative: ", x$alternative, "  estimated p value = ", x$est.p.value, " exact p value = ", x$exact.p.value,"\n")
  invisible(x)
}

# And a plotting function
plot.weighted_jt_statistic <- function(x, xlim=c(min(c(x$permutationen, x$statistic)), max(c(x$permutationen, x$statistic))), 
                                          col=c("gray","red"), main="JT permutation test",xlab="Statistic", ylab="Density", ...)
{
  plot(density(x$permutationen), col=col[1], xlim=xlim, main=main, xlab=xlab, ylab=ylab, ...)
  abline(v=x$statistic, col=col[2], lwd=3)
  legend("topright", col=col, legend=c("Permutations","Reference"))
  invisible(x)
}