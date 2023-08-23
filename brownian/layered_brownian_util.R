consc = function(from, to, by = 1)
{
  if ((to-from)*by < 0)
  {
    em = vector()
    return(em)
  }
  return(seq(from,to,by=by))
}

sigma_term = function(s,u1,u2,K,j)
{
  t1 = exp(-2/s*(2*K*j - K - u1) * (2*K*j - K - u2))
  t2 = exp(-2/s*(2*K*j - K + u1) * (2*K*j - K + u2))
  return(t1+t2)
}

tau_term = function(s,u1,u2,K,j)
{
  t1 = exp(-2*j/s* (4*K^2*j + 2*K*(u1-u2)))
  t2 = exp(-2*j/s* (4*K^2*j + 2*K*(u2-u1)))
  return(t1+t2)
}

check_gamma = function(u,total_time,start,end,limit)
{
  if (max(abs(start),abs(end)) > limit)
  {
    return(FALSE)
  }
  j = 1
  lb = 1 - sigma_term(total_time, start, end, limit, j)
  ub = lb + tau_term(total_time, start, end, limit, j)
  while (TRUE)
  {
    if (u < lb)
      return(TRUE)
    if (u > ub)
      return(FALSE)
    j = j+1
    lb = lb - sigma_term(total_time, start, end, limit,j)
    ub = lb + tau_term(total_time, start, end, limit,j)
  }
}

fetch_bound = function(n, ai, d)
{
  if (n < 1)
  {
    return(0)
  }
  l = length(ai)
  if (n <= l)
  {
    return(ai[n])
  } else
  {
    return(ai[l]+d*(n-l))
  }
}

ksi_term = function(s,u2,K,j)
{
  return((2*K*j+u2)*exp(-2*K*j*(K*j+u2)/s))
}

check_delta = function(total_time, start, end, bound1, bound2 = FALSE)
{
  if (max(abs(start),abs(end))> bound1)
  {
    return(FALSE)
  }
  if (end == 0)
  {
    end = start; start = 0
  }
  if (start == 0)
  {
    if (bound2 == FALSE)
    {
      return(check_delta_zero_1(total_time, end, bound1))
    } else
    {
      return(check_delta_zero_2(total_time, end, bound1, bound2))
    }
  } else
  {
    if (bound2 == FALSE)
    {
      return(check_delta_1(total_time, start, end, bound1))
    } else
    {
      return(check_delta_2(total_time, start, end, bound1, bound2))
    }
  }
}

check_delta_zero_1 = function(total_time, end, bound1)
{
  u = runif(1)
  count = 1
  lb = 1 - ksi_term(total_time, -end, bound1, count)/end
  ub = lb + ksi_term(total_time, end, bound1, count)/end
  while(TRUE)
  {
    if (u < lb)
      return(TRUE)
    if (u > ub)
      return(FALSE)
    count = count + 1
    lb = lb - ksi_term(total_time, -end, bound1, count)/end
    ub = lb + ksi_term(total_time, end, bound1, count)/end
  }
}

check_delta_zero_2 = function(total_time, end, bound1, bound2)
{
  u = runif(1)
  count = 1
  lb_n = end - ksi_term(total_time, -end, bound1, count)/end
  lb_d = end - ksi_term(total_time, -end, bound2, count)/end
  ub_n = lb_n + ksi_term(total_time, end, bound1, count)/end
  ub_d = lb_d + ksi_term(total_time, end, bound2, count)/end
  while(TRUE)
  {
    if (u < lb_n/ub_d)
      return(TRUE)
    if (u > ub_n/lb_d)
      return(FALSE)
    count = count + 1
    lb_n = lb_n - ksi_term(total_time, -end, bound1, count)/end
    lb_d = lb_d - ksi_term(total_time, -end, bound2, count)/end
    ub_n = lb_n + ksi_term(total_time, end, bound1, count)/end
    ub_d = lb_d + ksi_term(total_time, end, bound2, count)/end
  }
}

check_delta_1 = function(total_time, start, end, bound1)
{
  u = runif(1)
  j = 1
  Z = 1/(1-exp(-2*start*end/total_time))
  limit = bound1*0.5
  start = start - limit; end = end - limit
  lb = (1 - sigma_term(total_time, start, end, limit, j))*Z
  ub = lb + tau_term(total_time, start, end, limit, j)*Z
  while (TRUE)
  {
    if (u < lb)
      return(TRUE)
    if (u > ub)
      return(FALSE)
    j = j+1
    lb = lb - sigma_term(total_time, start, end, limit,j)*Z
    ub = lb + tau_term(total_time, start, end, limit,j)*Z
  }
}

check_delta_2 = function(total_time, start, end, bound1, bound2)
{
  u = runif(1)
  j = 1
  limit1 = bound1*0.5
  limit2 = bound2*0.5
  u1_K = start - limit1; u2_K = end - limit1
  u1_L = start - limit2; u2_L = end - limit2
  lb_n = 1 - sigma_term(total_time, u1_K, u2_K, limit1, j)
  lb_d = 1 - sigma_term(total_time, u1_L, u2_L, limit2, j)
  ub_n = lb_n + tau_term(total_time, u1_K, u2_K, limit1, j)
  ub_d = lb_d + tau_term(total_time, u1_L, u2_L, limit2, j)
  while (TRUE)
  {
    if (u < lb_n/ub_d)
      return(TRUE)
    if (u > ub_n/lb_d)
      return(FALSE)
    j = j+1
    lb_n = lb_n - sigma_term(total_time, u1_K, u2_K, limit1, j)
    lb_d = lb_d - sigma_term(total_time, u1_L, u2_L, limit2, j)
    ub_n = lb_n + tau_term(total_time, u1_K, u2_K, limit1, j)
    ub_d = lb_d + tau_term(total_time, u1_L, u2_L, limit2, j)
  }
}
