# working script for gamma diversity from dropbox
# just to get a handle on the data


shallow_total <- shallow %>%
  select(-Group)  %>%
  summarise_all(sum) 

shallow_total[shallow_total==0] <- NA
shallow_total <- shallow_total[!is.na(shallow_total)]


deep_total <- deep %>%
  select(-Group,-`Sub-Group`, -Site) %>% #drops the columns group and sub.group
  summarise_all(sum)

deep_total[deep_total==0] <- NA
deep_total <- deep_total[!is.na(deep_total)]



#basic stats
#look at the diversity on the Expected Number of Species at 400 individuals!!!!
#it does the flip flop you thought could occur

gamma_div <- data.frame(diversity(shallow_total))
colnames(gamma_div )[1] <- "H"
gamma_div$simp <- diversity(shallow_total, "simpson")
gamma_div$rich <- specnumber(shallow_total)
gamma_div$abund <-sum(shallow_total)
gamma_div$Group <- "Shallow"
gamma_div $rare400 <- rarefy(shallow_total, 400)


gamma_div[2,'H'] <-diversity(deep_total)
gamma_div[2,'simp'] <- diversity(deep_total, "simpson")
gamma_div[2,'rich'] <- specnumber(deep_total)
gamma_div[2,'abund'] <-sum(deep_total)
gamma_div[2,'Group'] <- "Deep"
gamma_div[2,'rare400'] <- rarefy(deep_total, 400)


#rarefaction


#rarefaction analysis and plot
out1 <- iNEXT(deep_total, q=0, datatype="abundance")
ggiNEXT(out1, type=1) + theme_bw(base_size=10) + xlim(0,700) + ylim (0,80)

out2 <- iNEXT(shallow_total, q=0, datatype="abundance")
ggiNEXT(out2, type=1) + theme_bw(base_size=10) + xlim(0,700) + ylim (0,80)


#get everything on the same plot with even more tedious and ridicolous code
#this is not worth explaining other than it extracts out the bits of data needed from the
#mess that iNext spits out so we can plot it how would like.

df_deep <- fortify(out1, type=1)
df_shallow <- fortify(out2, type=1)

df.point_deep <- df_deep[which(df_deep$method=="observed"),]
df.line_deep <- df_deep[which(df_deep$method!="observed"),]
df.line_deep$method <- factor(df.line_deep$method, 
                              c("interpolated", "extrapolated"),
                              c("interpolation", "extrapolation"))

df.point_shallow <- df_shallow[which(df_shallow$method=="observed"),]
df.line_shallow <- df_shallow[which(df_shallow$method!="observed"),]
df.line_shallow$method <- factor(df.line_shallow$method, 
                                 c("interpolated", "extrapolated"),
                                 c("interpolation", "extrapolation"))



ggplot(df_deep, aes(x=x, y=y)) + 
  ggtitle("Diversity of Systems (Gamma)") +
  geom_point(aes(shape=site), size=5, data=df.point_deep, color="blue3") + #plots the end points of the curves
  geom_line(aes(linetype=method), lwd=1.5, data=df.line_deep, color="blue3") + #plots acutal lines
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,
                  fill=site, colour=NULL), alpha=0.2, fill="blue3") + #plots envelopes
  #repeating the same for shallow
  geom_point(aes(shape=site), size=5, data=df.point_shallow, color="dodgerblue") +
  geom_line(aes(linetype=method), lwd=1.5, data=df.line_shallow, color="dodgerblue") +
  geom_ribbon(data=df_shallow, aes(ymin=y.lwr, ymax=y.upr,
                                   fill=site, colour=NULL), alpha=0.2, fill="dodgerblue") +
  labs(x="Number of Individuals", y="Species diversity") + #labelling the axis
  theme_bw(base_size=14) +
  theme(legend.position="right")
