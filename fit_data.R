dataframe = read.csv('ReedWordScore.csv', row.names = 1)

# fit gamma to each word frequency
column_names = c('strong_masculine', 'weak_masculine', 'strong_feminine', 'weak_feminine')
shapes = c()
rates = c()
for (i in 1:4)
{
  cmean = mean(dataframe[[column_names[i]]])
  cvar = var(dataframe[[column_names[i]]])
  rates[i] = cmean/cvar
  shapes[i] = cmean^2/cvar
}

print(rates)
print(shapes)

# aggregate scores
lambda =2
dataframe$score = 1*dataframe$Gender_Score_Gaucher2 + dataframe$Embedding_score*lambda

# linear regression on words
lm_result = lm(score ~ strong_masculine + weak_masculine + strong_feminine + weak_feminine, data = dataframe)
print(summary(lm_result))
coefs = lm_result$coefficients

intercept = coefs[1]
beta = coefs[2:5]
print(cor(cbind(dataframe$score, t(beta%*%t(dataframe[,2:5]))+intercept)))