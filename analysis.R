library(ggplot2)
setwd("~/Desktop/Eisai/Eisai/")
df = read.csv('table.csv')

ggplot(df[df$Group=='Jensen_Eisai_groupA',], aes(x=Duration,color=Group)) + geom_histogram()
ggplot(df[df$Layer=='LL-global_detected_seizures',], aes(x=Channels_1)) + geom_histogram() + facet_grid(Channels_1 ~Group)


tmpdatavg = ddply(tmpdat, .(subjid,time),summarize,meanfeat=mean(feat,na.rm=TRUE),sd=sd(feat,na.rm=TRUE),cofvar=sd/meanfeat)
#ggplot(tmpdatavg,aes(x=time,y=sd/meanfeat)) + geom_point() + geom_smooth(method='lm') + ylab(paste('log(',featName,')')) + xlab('Time (Day)') + ggtitle(paste('Individual fits for ',featName,' by subject')) + scale_x_continuous(labels = formatter2)  + facet_wrap(~subjid,scales="free")

g = ggplot(tmp[(tmp$feature %in% powerFeats),],aes(x=time,y=meanfeat)) + stat_smooth(aes(colour="Individual Trends",group=subjid),method="loess",span=smoothSpan,se=F,size=0.4,alpha=0.6,level=0.95)    + stat_summary(fun.data=mean_se,geom="ribbon",alpha=0.35)    + stat_smooth(aes(colour="Mean"),se=F,alpha=0.8,level=0.95)  + facet_wrap(~feature) +
  geom_hline(aes(colour="Reference Threshold",yintercept=upper_thres),size=1,linetype="dotted")  + theme(axis.text = element_text(size=12),axis.title = element_text(size=14),strip.text.x = element_text(size = 14),legend.position='top',legend.text=element_text(size=12)) + coord_cartesian(ylim = c(0.2,1),xlim=c(0,MAXTIMEPLOT)) +
  xlab('Time (Days)') + ylab('Normalized Feature Value') + ggtitle('Mean value - Spectral Features') + scale_colour_manual(name="",values=c("Individual Trends"="orange","Mean" = "blue","Reference Threshold"="red"))
g = g + guides(colour = guide_legend(override.aes = list(linetype = c(1,1,3))))
ggsave(paste('allpower_meanfeat',MAXTIMEPLOT,'.pdf',sep=''),dpi=600,height=6,width=8)

