# The statistic is calculated under the precondition that for each individual id the value x is known but the group g is not clear
# therefore. This is important for 
# It is not checked whether 

weighted.jonckheere.test <- function(x, g, wt, id = factor(x), B=1e4, normalize.wts = FALSE, alternative = c("two.sided", "increasing", "decreasing"),
                                     grouping = c("none","by.groups","by.values"))
{
  # check argument
  alternative = match.arg(alternative, c("two.sided", "increasing", "decreasing"))
  grouping = match.arg(grouping, c("none","by.groups","by.values"))

  # Clean factor "id"
  id <- factor(as.character(id))
  
  # Normalize
  if (normalize.wts) {
    wt = wt / tapply(wt, id, sum)[id]
  }
  # First, we calculate the statistic under the given condition
  basis = .Call("rcpp_wt_jonckheere", values = as.double(x), groups = as.integer(g), weights = as.double(wt), PACKAGE="wtJonckheere")
    
   
  # Select between sampling types
  if (grouping == "by.groups") 
  {  
    # do a check: x must be identical for all groups
    check.v = tapply(x, id, max)
    if (any( x != check.v[id])) stop("Error: x values must be identical for each id individual!\n");
    
    # calculate multiple permutations. For this purpose, the values for x are exchanged between the ids
    perm_x_values = tapply(x, id, mean)
    permutationen <- sapply(1:B, function(dum) {
      perm_x_new_values = perm_x_values[order(runif(length(perm_x_values)))]
      names(perm_x_new_values) <- names(perm_x_values)
      return(.Call("rcpp_wt_jonckheere", values = as.double(perm_x_new_values[id]), groups = as.integer(g), weights = as.double(wt), PACKAGE="wtJonckheere"))
    })
  }  
  if (grouping == "by.values") 
  {  
    # do a check: g must be identical for all groups
    check.v = tapply(g, id, max)
    if (any( g != check.v[id])) stop("Error: group values must be identical for each id individual!\n");
    
    # calculate multiple permutations. For this purpose, the values for x are exchanged between the ids
    perm_g_values = tapply(g, id, mean)
    permutationen <- sapply(1:B, function(dum) {
      perm_g_new_values = perm_g_values[order(runif(length(perm_g_values)))]
      names(perm_g_new_values) <- names(perm_g_values)
      return(.Call("rcpp_wt_jonckheere", values = as.double(x), groups = as.integer(perm_g_new_values[id]), weights = as.double(wt), PACKAGE="wtJonckheere"))
    })
  }  
  if (grouping =="none")
  {
    permutationen <- sapply(1:B, function(dum) {
      return(.Call("rcpp_wt_jonckheere", values = as.double(x[order(runif(length(x)))]), 
                                         groups = as.integer(g[order(runif(length(g)))]), weights = as.double(wt), PACKAGE="wtJonckheere"))
    })
  }
    
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
  
  result = list(statistic = basis, permutations=permutationen, alternative = alternative, 
                mean.permutation = mittelwert, std.permutation = std, n.permutation=B,
                est.p.value = pv, exact.p.value = exakt.pv);
  class(result) = "weighted_jt_statistic"
  return(result)
}

# And a printing function
summary.weighted_jt_statistic <- function(x, ...)
{
  cat("\nWeighted Jonckheere Terpstra Testing\n\n")
  cat("Alternative: ", x$alternative, "\n")
  cat("JT statistic of data = ", x$statistic,"\n")
  cat("Summary of ", x$n.permutation, " permutations. Average JT statistic =",signif(x$mean.permutation,3), " +/- std. dev. ", signif(x$std.permutation,3), "\n")
  qq = signif(quantile(x$permutations, probs=c(0.025,0.5,0.975)),3)
  cat("2.5% quantile = ", qq[1], " median= ", qq[2], " 97.5% quantile = " ,qq[3],"\n\n")
  cat("estimated p value = ", signif(x$est.p.value,5), " \t exact p value = ", signif(x$exact.p.value,5),"\n")
  invisible(x)
}

# And a plotting function
plot.weighted_jt_statistic <- function(x, xlim=c(min(c(x$permutations, x$statistic)), max(c(x$permutations, x$statistic))), 
                                          col=c("gray","red",rgb(1,0,0,0.5)), main="JT permutation test",xlab="Statistic", ylab="Frequency", ...)
{
 
  breite = max(xlim)-min(xlim)
#   dens <- density(x$permutations);
  histx = hist(x$permutations, plot=F)
  faktor = mean(diff(histx$breaks))*x$n.permutation
  
  mitte = max(histx$counts)/3
  
  plot(histx, col=col[1], xlim=xlim, main=main, xlab=xlab, ylab=ylab, ...)
  xe <- seq(min(xlim), max(xlim), len=1000)
  lines(faktor * dnorm(xe, mean=x$mean.permutation, sd=x$std.permutation) ~ xe, col=col[3], lwd=3);
  abline(v=x$statistic, col=col[2], lwd=3)
  text(x=median(x$permutations), y=mitte, label="Permutations", col=col[1])
  text(x=median(x$permutations), y=mitte*1.1, label="Std. normal approximation", col=col[3])
  text(x=x$statistic - breite/20, y=mitte, col=col[2], label="Reference", srt=90)
  invisible(x) 
}