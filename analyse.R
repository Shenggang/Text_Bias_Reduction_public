library(ggplot2)
library(dplyr)
library(gridExtra)

bias_orig = as.numeric(dataframe$score)
bias_before = as.numeric(beta%*%t(dataframe[,2:5])) + intercept
bias_after = as.numeric(beta%*%t(round(result.mean))) + intercept


corr = cor(cbind(bias_orig, bias_before))
cat("Correlation between orignal score and fitted score is", corr[1,2], '\n')

improv = abs(bias_before)-abs(bias_after)
cat('mean improvement =', mean(improv), '\n')

df_all = data.frame(before=bias_before, after=bias_after, improv=improv, positive=(improv>=0))
df_pos = df_all[df_all$positive, ]
per = 0.75
df_pos$improv75 = (df_pos$improv/abs(df_pos$before))>0.75
df_65 = df_all[abs(df_all$before)>per, ]
df_20 = df_all[abs(df_all$before)>0.23, ]

#width 900 height 300
cat('mean improvement percentage =', mean(improv/abs(bias_before)*100), '% \n')
#cat('mean positive improvement =', mean(pos_res[,1]), '\n')

#cat('mean positive improvement percentage =', mean(pos_res[,1]/abs(pos_res[,2])*100), '% \n')

#cat('mean abs negative bias before = ', mean(abs(neg_res[,2])),'\n')
cat('max negative before = ', max(abs(df_all[!df_all$positive, ]$before)), '\n')
#cat('mean abs bias before =', mean(abs(neg_res[,2])),'\n')

cat('positive percentage = ', length(df_pos$before)/length(df_all$before), '\n')
cat('75 improv percentage = ', sum(df_pos$improv75)/length(df_pos$before), '\n')
cat("mean before = ", paste(mean(abs(df_all$before)), mean(abs(df_pos$before)), 
                            mean(abs(df_20$before)), mean(abs(df_65$before)),  sep=' & '), '\n')

cat("mean after = ", paste(mean(abs(df_all$after)), mean(abs(df_pos$after)), 
                           mean(abs(df_20$after)) , mean(abs(df_65$after)), sep=' & '), '\n')

cat("mean improv = ", paste(mean(df_all$improv), mean(df_pos$improv), 
                            mean(df_20$improv), mean(df_65$improv),  sep=' & '), '\n')
cat("mean % improv = ", paste(mean(df_all$improv/abs(df_all$before)), mean(df_pos$improv/abs(df_pos$before)),
                              mean(df_20$improv/abs(df_20$before)), mean(df_65$improv/abs(df_65$before)), sep=' & '), '\n')

set.seed(1050)
# # plot randomly sampled points
df_sampled = sample_n(df_all, 3000)
plot1 = ggplot(df_sampled, aes(x=abs(before), y=improv, group=positive)) + geom_point(aes(color=positive, shape=positive)) +
  scale_color_manual(values=c('#EE0000','#00FAFA')) + scale_shape_manual(values=c(16, 3)) + scale_size_manual(values=c(1,2))+
  labs(title="(A) Improvement vs. bias before mitigation", x="Bias before mitigation", y="Raw improvement score") +
  geom_abline(intercept = 0, slope=1, color="darkgrey", linetype="dashed", size=0.9) +
  theme(text=element_text(size=15))
# plot percentage improvment positive only
df_pos_sampled = sample_n(df_pos, 3000)
plot2 = ggplot(df_pos_sampled, aes(x=abs(before), y=improv/abs(before)*100, group=improv75)) + 
  geom_point(aes(color=improv75, shape=improv75)) +
  scale_color_manual(values=c('#EE0000','#00FAFA')) + scale_shape_manual(values=c(16, 3)) + scale_size_manual(values=c(1,2))+
  labs(title="(B) Percentage (positive) improvement", x="Bias before mitigation", y="Improvement percentage")+
  theme(text=element_text(size=15))

width=12
height=8
g1 = ggplotGrob(plot1)
g2 = ggplotGrob(plot2)
grid::grid.newpage()
grid::grid.draw(rbind(g1,g2))
g = arrangeGrob(plot1, plot2, nrow=2)
#ggsave('improv_subfigures.eps', g, width=width, height=height, device='eps', dpi=300)
ggsave('improv_subfigures.jpg', g, width=width, height=height, device='jpg', dpi=300)

# plot histogram of bias before and after
hist1 = ggplot(df_all, aes(x=before))+
      geom_histogram(color='darkgreen', fill='lightgreen', binwidth=0.1)+
      scale_x_continuous(breaks=seq(-4,4,0.5), lim=c(-4,4))+
      labs(title="(A) Bias before mitigation", x="Bias before mitigation", y="Count") +
        theme(text=element_text(size=15))
hist2 = ggplot(df_all, aes(x=after))+
      geom_histogram(color='darkblue', fill='darkred',binwidth = 0.05)+
      scale_x_continuous(breaks=seq(-4,4,0.5), lim=c(-4,4))+
      labs(title="(B) Bias after mitigation", x="Bias after mitigation", y="Count") +
      theme(text=element_text(size=15))
g1 = ggplotGrob(hist1)
g2 = ggplotGrob(hist2)
grid::grid.newpage()
grid::grid.draw(rbind(g1,g2))
g = arrangeGrob(hist1, hist2, nrow=2)
#ggsave('bias_plots.eps', width=width, height=height, device='eps',dpi=300)
ggsave('bias_plots.jpg', g, width=width, height=height, device='jpg', dpi=300)

