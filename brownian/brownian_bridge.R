brownian_bridge = function(total_time, start, end, time_points)
{
  if (length(time_points) == 0)
  {
    return(rbind(c(0,total_time),c(start, end)))
  }
  time_points = sort(time_points)
  if (tail(time_points, n=1) >= total_time)
  {
    print("Required time point exceeds total time. Unable to simulate Brownian bridge")
    return()
  }
  bm = brownian_motion(total_time, time_points)
  bm_sq = bm[2,]
  time_points = c(0, time_points, total_time)
  r_time = time_points/total_time
  bb = bm_sq - r_time*bm_sq[length(bm_sq)]
  return(rbind(bm[1,],start*(1-r_time) + end * r_time + bb))
}

brownian_bridge_with_min = function(total_time, start, end, min, min_t, time_points)
{
  if (length(time_points) == 0)
  {
    return(rbind(c(0, min_t, total_time),c(start, min, end)))
  }
  time_points = sort(time_points)
  if (tail(time_points, n=1) >= total_time)
  {
    print("Required time point exceeds total time. Unable to simulate Brownian bridge")
    return()
  }
  th1 = sum(time_points<min_t)
  th2 = sum(time_points<=min_t)+1
  bb1_xt = time_points[consc(1,th1)]
  bb2_xt = time_points[consc(th2,length(time_points))]
  beb1 = bessel_bridge(min_t, start - min, min_t-rev(bb1_xt))
  beb2 = bessel_bridge(total_time - min_t, end - min, bb2_xt-min_t)
  bb = c(rev(beb1[2,])+min, beb2[2,2:length(beb2[2,])]+min)
  ts = c(0,bb1_xt, min_t, bb2_xt, total_time)
  return(rbind(ts,bb))
}

brownian_bridge_with_max = function(total_time, start, end, max, max_t, time_points)
{
  bb = brownian_bridge_with_min(total_time, -start, -end, -max, max_t, time_points)
  bb[2,] = -bb[2,]
  return(bb)
}

r_min_brownian_bridge = function(total_time, start, end, interval = FALSE)
{
  y = end - start
  if (length(interval) == 1)
  {
    E = rexp(1)
  } else
  {
    interval = sort(interval) - start
    d = ((y-2*interval[1])^2-y^2)/2/total_time
    c = ((y-2*interval[2])^2-y^2)/2/total_time
    E = rexp_interval(1,c,d)
  }
  Z1 = 0.5*(y - sqrt(2*total_time*E+y^2))
  c1 = (y-Z1)^2/2/total_time
  c2 = Z1^2/2/total_time
  I1 = rinvgauss(1, mean = sqrt(c1/c2), dispersion = 2*c1)
  I2 = 1/rinvgauss(1, mean = sqrt(c2/c1), dispersion = 2*c2)
  if (runif(1) < 1/(1+sqrt(c1/c2)))
  {
    V = I1
  } else
  {
    V = I2
  }
  Z2 = total_time/(1+V)
  return(c(Z1+start, Z2))
}

r_max_brownian_bridge = function(total_time, start, end, interval=FALSE)
{
  result = r_min_brownian_bridge(total_time, -start, -end, -interval)
  return(c(-result[1], result[2]))
}
