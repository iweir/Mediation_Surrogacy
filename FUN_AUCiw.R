AUCiw <- function(x, y, na.rm = FALSE, to=max(x), from=min(x)) {
  
  if (length(x) != length(y))
    stop("length x must equal length y")
  
  values   <- approx(x, y, xout = sort(unique(c(from, to, x[x > from & x < to]))))
  idx      <- order(values$x)
  values$x <- values$x[idx]
  values$y <- values$y[idx]
  
  # use step function summation
  a        <- sum( values$y[-length(values$y)] * (values$x[-1] - values$x[-length(values$x)]))
  return(a)
}

## adapted from the AUC() function in the DescTools library and auc() function in the MESS library
