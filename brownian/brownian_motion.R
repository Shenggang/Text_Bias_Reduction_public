brownian_motion = function(total_time, time_points) 
{
  if (length(time_points) == 0)
  {
    return(rbind(c(0, total_time),c(0, rnorm(1)*sqrt(total_time))))
  }
  time_points = sort(time_points)
  if (tail(time_points, n=1) >= total_time)
  {
    print("Required time point exceeds total time. Unable to simulate Brownian motion")
    return()
  }
  time_points = c(0, time_points, total_time)
  time_diff = diff(time_points)
  z = rnorm(length(time_diff))*sqrt(time_diff)
  bm_sq = c(0, cumsum(z))
  return(rbind(time_points,bm_sq))
}
