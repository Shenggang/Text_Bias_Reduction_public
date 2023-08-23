total_entries = length(summary)

result.mean = matrix(0, nrow=total_entries, ncol=4)
result.ESS = c()

for (i in 1:total_entries)
{
  result.mean[i, ] = summary[[i]]$r.mean
  result.ESS[i] = summary[[i]]$r.ESS
}