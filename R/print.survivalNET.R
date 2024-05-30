
print.survivalNET <- function (x, digits=4, ...)
{
  print(round(x$t.table, digits = digits))
  
  lrs <- 2 * (x$loglik[1] - x$loglik[2])
  df <- length(gsub("\\+", "", attr(terms(x$formula), "term.labels")))
  pv <- 1 - pchisq(q=lrs, df=df)
  nmiss <- sum(x$missing)
  cat("\n", "Likelihood ratio test=", round(lrs, digits),
      " on ", df, " df, p=", round(pv, digits), sep="")
  
  cat("\n", "n=", x$n, ", ", "number of events=", x$nevent, sep="")
  
  if(nmiss==1) { cat("\n", "(", nmiss, " observation deleted due to missingness)", sep="") }
  if(nmiss >1) { cat("\n", "(", nmiss, " observations deleted due to missingness)", sep="") }
}



