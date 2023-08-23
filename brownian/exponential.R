rexp_interval = function(n,lowerbound, upperbound)
{
  result = rep(0,n)
  c = lowerbound; d = upperbound
  if (d-c <= 500)
  {
    for (i in 1:n)
    {
      y = runif(1)
      lz = d - log(exp(d-c)-1)
      r = lz - log(y)
      result[i] = r - log(exp(r-c)-1)
    }
  } else 
  {
    if (c > 700)
      print("Interval too large or lowerbound too large, likely not able to return value")
    for (i in 1:n)
    {
      y = runif(1)
      Z = 1/(exp(-c) - exp(-d))
      result[i] = -log(exp(-c) - y/Z)
    }
  }
  return(result)
}