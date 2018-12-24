setwd("~/Dropbox/R/transients/githubupload/rev3/")
setwd("C:/Users/Hoameng/Dropbox/R/transients/githubupload/rev3/")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
packages= c("numDeriv","RCurl","knitr","MASS","ggplot2","orcutt","reshape2","plyr",'RColorBrewer','R.matlab','zoo','data.table','scales',"lme4","nlme","car","lmtest","coefplot2","GGally","scales",'Rmisc',"tidyr","effsize","boot",'grid','geepack','gee')
lapply(packages,usePackage)
smoothSpan = 0.65

#### Remove Outlier Function ####
#Function to remove outliers from raw data. These outliers were manually identified due to large perturbations
#in features, that were subsequently validated to be artifactual as described in the manuscript.
removeOutliers = function(dat) {
  tmpdat = dat
  tmpdat$subjid = factor(tmpdat$subjid)
  tmpdat$ch = factor(tmpdat$ch)
  tmpdat$szCh = factor(tmpdat$szCh)
  # idx = tmpdat$subjid=="I004_A0001_D00"
  # 
  # idx = tmpdat$subjid=="I004_A0001_D00" & (tmpdat$time>42)
  # tmpdat[idx,]$feat = NA
  # 
  # idx = tmpdat$subjid=="I004_A0002_D00"
  # 
  # idx = tmpdat$subjid=="I004_A0002_D00" & (tmpdat$time>100 & tmpdat$time<250)
  # tmpdat[idx,]$feat = NA
  # 
  # #A3
  # idx = tmpdat$subjid=="I004_A0003_D00"
  # 
  # #
  # idx = tmpdat$subjid=="I004_A0004_D00"
  # 
  # idx = tmpdat$subjid=="I004_A0004_D00" & tmpdat$time>230
  # tmpdat[idx,]$feat = NA

  
  #23_002
  idx = tmpdat$subjid=="NVC1001_23_002" 
  
  idx = tmpdat$subjid=="NVC1001_23_002" & (tmpdat$ch==16 & tmpdat$time>100 & tmpdat$time < 200)
    tmpdat[idx,]$feat = NA
  
  idx = tmpdat$subjid=="NVC1001_23_002" & (tmpdat$time>20 & tmpdat$time < 60)
    tmpdat[idx,]$feat = NA
  
  ##23_003
  idx = tmpdat$subjid=="NVC1001_23_003" 

  idx = tmpdat$subjid=="NVC1001_23_003" & (tmpdat$time>485 & tmpdat$time < 530) & (tmpdat$ch %in% c(13,14,15,16))
  tmpdat[idx,]$feat = NA
  idx = tmpdat$subjid=="NVC1001_23_003" & (tmpdat$time>640 & tmpdat$time < 694)
  tmpdat[idx,]$feat = NA
  idx = tmpdat$subjid=="NVC1001_23_003" & (tmpdat$ch==8 & tmpdat$time > 50 & tmpdat$time < 75)
  tmpdat[idx,]$feat = NA
  
  #23_004
  idx = tmpdat$subjid=="NVC1001_23_004" 

  
  #23_005
  idx = tmpdat$subjid=="NVC1001_23_005" 

  idx = tmpdat$subjid=="NVC1001_23_005" & (tmpdat$time>215 | tmpdat$time<=1)
    tmpdat[idx,]$feat = NA
  idx = tmpdat$subjid=="NVC1001_23_005" & (tmpdat$time>125 & tmpdat$ch==1)
    tmpdat[idx,]$feat = NA
  idx = tmpdat$subjid=="NVC1001_23_005" & (tmpdat$time>170 & tmpdat$ch==8)
    tmpdat[idx,]$feat = NA
  idx = tmpdat$subjid=="NVC1001_23_005" & tmpdat$time==4
    tmpdat[idx,]$feat = NA

  
  #23_006
  idx = tmpdat$subjid=="NVC1001_23_006" 

  idx = tmpdat$subjid=="NVC1001_23_006" & (tmpdat$time<5)
  tmpdat[idx,]$feat = NA
  idx = tmpdat$subjid=="NVC1001_23_006" & (tmpdat$time>40 & tmpdat$time<43)
  tmpdat[idx,]$feat = NA
  idx = tmpdat$subjid=="NVC1001_23_006" & (tmpdat$time>=229 & tmpdat$time<=260) & (tmpdat$ch %in% c(5,6))
  tmpdat[idx,]$feat = NA
  
  #23_007   
  idx = tmpdat$subjid=="NVC1001_23_007" 

  idx = tmpdat$subjid=="NVC1001_23_007" & (tmpdat$time<2) #dropouts at beginning of recording
  tmpdat[idx,]$feat = NA
  
  #24_001
  idx = tmpdat$subjid=="NVC1001_24_001" 
  
  #24_002   
  idx = tmpdat$subjid=="NVC1001_24_002" 

  idx = tmpdat$subjid=="NVC1001_24_002" & (tmpdat$time>25 & tmpdat$time<50) & (tmpdat$ch %in% c(13,16))
  tmpdat[idx,]$feat = NA
  idx = tmpdat$subjid=="NVC1001_24_002" & (tmpdat$time>50 & tmpdat$time<80) & (tmpdat$ch %in% c(2,3))
  tmpdat[idx,]$feat = NA
  idx = tmpdat$subjid=="NVC1001_24_002" & (tmpdat$time==325 | tmpdat$time==382 | tmpdat$time == 495)
  if (sum(idx)>0) {
    tmpdat[idx,]$feat = NA
  }
  
  #24_004 
  idx = tmpdat$subjid=="NVC1001_24_004" 

  idx = tmpdat$subjid=="NVC1001_24_004" & (tmpdat$ch==10 & ((tmpdat$time>43 & tmpdat$time<=71) ))
  tmpdat[idx,]$feat = NA
  
  idx = tmpdat$subjid=="NVC1001_24_004" & (tmpdat$ch==1 & (tmpdat$time>140 & tmpdat$time<150))
  tmpdat[idx,]$feat = NA
  
  idx = tmpdat$subjid=="NVC1001_24_004" & (tmpdat$ch==4 & (tmpdat$time>75 & tmpdat$time<153))
  tmpdat[idx,]$feat = NA
  
  idx = tmpdat$subjid=="NVC1001_24_004" & (tmpdat$time>125 & tmpdat$time<190) & (tmpdat$ch %in% c(10,11,12))
  tmpdat[idx,]$feat = NA
  
  idx = tmpdat$subjid=="NVC1001_24_004" & (tmpdat$time>295 & tmpdat$time<323) & (tmpdat$ch %in% c(1,2,4,15,16))
  tmpdat[idx,]$feat = NA
  
  idx = tmpdat$subjid=="NVC1001_24_004" & (tmpdat$time>325)
  tmpdat[idx,]$feat = NA
  
  idx = tmpdat$subjid=="NVC1001_24_004" & (tmpdat$time<2)
  tmpdat[idx,]$feat = NA
  #24_005 
  idx = tmpdat$subjid=="NVC1001_24_005" 
  
  #25_001 #a lot of variability from dropped signals
  idx = tmpdat$subjid=="NVC1001_25_001" 
  ggplot(tmpdat[idx,],aes(x=time,y=feat)) + geom_point() + facet_wrap(~ch)
  ggplot(tmpdat[idx,],aes(x=time,y=feat,colour=ch)) + geom_point() 
  idx = tmpdat$subjid=="NVC1001_25_001" & (tmpdat$time>18 & tmpdat$time<52) & (tmpdat$ch==7 | tmpdat$ch==12)
  tmpdat[idx,]$feat = NA
  
  #25_002
  idx = tmpdat$subjid=="NVC1001_25_002" 

  # idx = tmpdat$subjid=="NVC1001_25_002" & ((tmpdat$ch==2 | tmpdat$ch==3 | tmpdat$ch==4 | tmpdat$ch==5 | tmpdat$ch == 6 | tmpdat$ch==7 | tmpdat$ch==8) & tmpdat$time>000)
  # tmpdat[idx,]$feat = NA
  idx = tmpdat$subjid=="NVC1001_25_002" & (tmpdat$time>49 & tmpdat$time<60)
  tmpdat[idx,]$feat = NA
  idx = tmpdat$subjid=="NVC1001_25_002" & (tmpdat$time>99 & tmpdat$time<141)
  tmpdat[idx,]$feat = NA
  idx = tmpdat$subjid=="NVC1001_25_002" & (tmpdat$time>440 & tmpdat$time<450)
  tmpdat[idx,]$feat = NA
  idx = tmpdat$subjid=="NVC1001_25_002" & (tmpdat$time %in% c(50, 55, 261, 390, 447))
  tmpdat[idx,]$feat = NA
  
  #25_003
  idx = tmpdat$subjid=="NVC1001_25_003" 

  
  #25_004
  idx = tmpdat$subjid=="NVC1001_25_004" 

  idx = tmpdat$subjid=="NVC1001_25_004" & (tmpdat$time<2)
  tmpdat[idx,]$feat = NA
  
  idx = tmpdat$subjid=="NVC1001_25_004" & (tmpdat$time>245 & tmpdat$time<350)
  tmpdat[idx,]$feat = NA

  idx = tmpdat$subjid=="NVC1001_25_004" & (tmpdat$time>390 & tmpdat$time<440)
  tmpdat[idx,]$feat = NA

  idx = tmpdat$subjid=="NVC1001_25_004" & (tmpdat$time>460 & tmpdat$time<550)
  tmpdat[idx,]$feat = NA
  # 
  #25_005
  idx = tmpdat$subjid=="NVC1001_25_005" 
  
  dogidx = (tmpdat$subjid == 'I004_A0001_D00' | tmpdat$subjid == 'I004_A0002_D00' | tmpdat$subjid == 'I004_A0003_D00' | tmpdat$subjid == 'I004_A0004_D00' )
  tmpdat = tmpdat[!dogidx,]
  
  return(tmpdat)
}

#### Utility Functions ####
summarySE = function(data=NULL, measurevar, groupvars=NULL, na.rm=TRUE, conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm),
                     median = median(xx[[col]],na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
normalizeCh = function(dat) {
  #normalize 95% each subject's chs to values of [0-1],
  subjects = unique(dat$subjid)
  for (i in 1:length(subjects)) {
    subj = subjects[i]
    feat = dat[dat$subjid==subj,]$feat
    ch = dat[dat$subjid==subj,]$ch
    uniqueCh = unique(ch)
    for (j in 1:length(uniqueCh)) {
      tmpfeat = feat[ch==uniqueCh[j]]
      q = quantile(tmpfeat,c(0.025,0.975),na.rm=TRUE)
      dat[dat$subjid==subj & dat$ch==uniqueCh[j],]$feat = (tmpfeat-q[1])/(q[2]-q[1])
    }
  }
  return(dat)
}
formatter2 = function(x){ 
  x
}

############# COF VAR ############
indivFit= function(dat,maxTime,featName) {
  tmp = dat[dat$time<maxTime,]
  tmpavg = ddply(tmp, .(subjid,time),summarize,meanfeat=mean(feat,na.rm=TRUE),sd=sd(feat,na.rm=TRUE),cofvar=sd/meanfeat)
  tmpavg = groupedData(meanfeat~time|subjid,data=tmpavg)
  mod = lmList(cofvar~time, data=tmpavg,na.action=na.omit)
  f1.lis = summary(mod)
  a = f1.lis$coefficients
  a = a[,,2]
  coeff = a[order(rownames(a)),]
  out = data.frame(subjid=rownames(coeff),feat=featName,b=coeff[,1],t=coeff[,3],p=coeff[,4])
  return(out)
}
countall = function(tmp) {
  noc = sum(tmp[,5]>0.05);
  pos = sum(tmp[tmp[,5]<=0.05,3]>0)
  neg = sum(tmp[tmp[,5]<=0.05,3]<0)
  return(c(noc,pos,neg))
}
cleanAndIndivFit = function(dat,featName,maxTime) {
  dat$subjid = factor(dat$subjid)
  dat$ch = factor(dat$ch)
  dat$szCh = factor(dat$szCh)
  str(dat)
  tmpdat = dat;
  tmpdat = removeOutliers(tmpdat)
  indivFits = indivFit(tmpdat,maxTime,featName)
  
  #plot fit
  tmpdatavg = ddply(tmpdat, .(subjid,time),summarize,meanfeat=mean(feat,na.rm=TRUE),sd=sd(feat,na.rm=TRUE),cofvar=sd/meanfeat)
  #ggplot(tmpdatavg,aes(x=time,y=sd/meanfeat)) + geom_point() + geom_smooth(method='lm') + ylab(paste('log(',featName,')')) + xlab('Time (Day)') + ggtitle(paste('Individual fits for ',featName,' by subject')) + scale_x_continuous(labels = formatter2)  + facet_wrap(~subjid,scales="free")
  #ggsave(paste('indiv_fits_cofvar',featName,'-',maxTime,'.pdf',sep=''),width=8,height=10,dpi=600)
  
  return(indivFits)
}

######### PLOT AVG FEATURES###########
indivFit= function(tmpdat,featName,maxTime,logflag) {
  tmp = tmpdat[tmpdat$time<maxTime,]
  tmpavg = ddply(tmp, .(subjid,time),summarize,feat=mean(feat))
  tmpavg = groupedData(feat~time|subjid,data=tmpavg)
  if (logflag ==0 ) {
    mod = lmList(feat~time, data=tmpavg,na.action=na.omit)
  }
  if (logflag ==1 ) {
    mod = lmList(log(feat)~time, data=tmpavg,na.action=na.omit)
  }
  if (logflag ==2 ) {
    mod = lmList(log(feat)~log(time), data=tmpavg,na.action=na.omit)
  }
  f1.lis = summary(mod)
  a = f1.lis$coefficients
  a = a[,,2]
  a = cbind(a,f1.lis$r.squared)
  colnames(a)[5] = 'r.squared'
  coeff = a[order(rownames(a)),]
  out = data.frame(subjid=rownames(coeff),feat=featName,b=coeff[,1],t=coeff[,3],p=coeff[,4],rsquared=coeff[,5])
  return(out)
}
indivFitAndPlot = function(tmpdat,featName,maxTime,logflag) {
  indivFits = indivFit(tmpdat,featName,maxTime,logflag)
  #plot fit
  tmpdatavg =  summarySE(tmpdat, measurevar="feat", groupvars=c("time","subjid"),na.rm=TRUE)
  if (logflag ==0 ) {
  #  ggplot(tmpdatavg,aes(x=time,y=feat)) + geom_point() + geom_smooth(method='lm') + ylab(featName) + xlab('Time (Day)') + ggtitle(paste('Individual fits for ',featName,' by subject')) + scale_x_continuous(labels = formatter2)  + facet_wrap(~subjid,scales="free")
  #  ggsave(paste('indiv_fits_nolog',featName,'-',maxTime,'.pdf',sep=''),width=8,height=10,dpi=600)  
  }
  if (logflag ==1 ) {
  #  ggplot(tmpdatavg,aes(x=time,y=log(feat))) + geom_point() + geom_smooth(method='lm') + ylab(paste('log(',featName,')')) + xlab('Time (Day)') + ggtitle(paste('Individual fits for ',featName,' by subject')) + scale_x_continuous(labels = formatter2)  + facet_wrap(~subjid,scales="free")
  #  ggsave(paste('indiv_fits_log',featName,'-',maxTime,'.pdf',sep=''),width=8,height=10,dpi=600)  
  }
  if (logflag ==2 ) {
  #  ggplot(tmpdatavg,aes(x=log(time),y=log(feat))) + geom_point() + geom_smooth(method='lm') + ylab(paste('log(',featName,')')) + xlab('log(Time (Day))') + ggtitle(paste('Individual fits for ',featName,' by subject')) + scale_x_continuous(labels = formatter2)  + facet_wrap(~subjid,scales="free")
  #  ggsave(paste('indiv_fits_loglog',featName,'-',maxTime,'.pdf',sep=''),width=8,height=10,dpi=600)
  }
  
  return(indivFits)
}
indivFits = NULL

id = '24hr2'
featShortList = c('LL','Area','HW','Energy','Delta','Theta','Alpha','Beta','Gamma','Gamma1','Gamma2','Higamma','higamma1','higamma2','firstder','secondder','RMS')
featLabelList = c('Line Length','Area','Halfwave Amp','Energy','Delta Power','Theta Power','Alpha Power','Beta Power','Gamma Power','Gamma (30-45 Hz)','Gamma (55-100 Hz)','High Gamma Power','High Gamma (100-140 Hz)','High Gamma (140-180 Hz)','First Derivative','Second Derivative','RMS')
timeFeats = featLabelList[c(1,2,3,4)]
derFeats = featLabelList[c(15,16,17)]
powerFeats = featLabelList[c(5,6,7,8,9,12)]
power2Feats = featLabelList[c(10,11,13,14)]

indivAndGroupPlot = function(featShort,featLabel,maxTime,logflag) {
  dat = as.data.frame(readMat(paste("combined_",featShort,"_",id,".mat",sep='')))
  #ggplot(dat[dat$subjid=='NVC1001_25_002',],aes(x=time,y=feat)) + geom_point() + facet_wrap(~ch,scales="free_y")
  #ggplot(dat[dat$subjid=='NVC1001_24_001',],aes(x=time,y=feat)) + geom_point() + facet_wrap(~ch,scales="free_y")
  
  #ggplot(dat[dat$subjid=='NVC1001_24_002',],aes(x=time,y=feat)) + geom_point() + facet_wrap(~ch,scales="free_y")
  
  #average channels

  #remove outliers
  dat = removeOutliers(dat)
  tmpdat_ch = summarySE(dat, measurevar="feat", groupvars=c("time","subjid"))
  tmpdat_ch$subjid = revalue(tmpdat_ch$subjid,c("NVC1001_23_002" = "1","NVC1001_23_003" = "2","NVC1001_23_004" = "3","NVC1001_23_005" = "4","NVC1001_23_006" = "5",
                                                      "NVC1001_23_007" = "6","NVC1001_24_001" = "7","NVC1001_24_002" = "8","NVC1001_24_004" = "9","NVC1001_24_005" = "10",
                                                      "NVC1001_25_001" = "11","NVC1001_25_002" = "12","NVC1001_25_003" = "13","NVC1001_25_004" = "14","NVC1001_25_005" = "15"))
  
  #ggplot(tmpdat_ch[tmpdat_ch$time<500,],aes(x=time,y=feat)) + geom_point() + facet_wrap(~subjid,scales="free_y")
  #ggsave(paste('subjectlevel-',featShort,500,'.pdf',sep=''),dpi=600,width=10,height=8)
  #ggplot(tmpdat_ch[tmpdat_ch$time<200,],aes(x=time,y=feat)) + geom_point() + facet_wrap(~subjid,scales="free_y")
  #ggsave(paste('subjectlevel-',featShort,200,'.pdf',sep=''),dpi=600,width=10,height=8)
  
  dat = dat[,-3]
  origdat = dat
  dat = dat[dat$time<=maxTime,]
  textsize = 18
  #before normalization, calculate Coefficient of variation, since this needs to be on a ratio scale and nonnegative.
  tmpdat_ch = summarySE(dat, measurevar="feat", groupvars=c("time","subjid"))
  
  #plot Coeff of Var between channels, averaged across all patients 
  tmpdat_ch$cofvar = tmpdat_ch$sd/tmpdat_ch$feat 
  tmpdat_meancofvarch = summarySE(tmpdat_ch, measurevar="cofvar", groupvars=c("time"))
  xlabel='Time (Day)'
  ylabel=paste('Mean CoV - ',featLabel,sep='')
  title=paste('Mean CoV across channels, averaged across subjects',sep='')
  g1 = ggplot(dat=tmpdat_meancofvarch, aes(x=time, y=cofvar)) + 
    geom_errorbar(aes(ymin=cofvar-se, ymax=cofvar+se), width=.1) +
    theme(legend.title=element_blank(),text = element_text(size=textsize),axis.text.x = element_text(angle=90, vjust=1)) + xlab(xlabel) + ylab(ylabel) + scale_x_continuous(labels = formatter2) 
  g1
  ggsave(paste('MeanCoV-subj-',featShort,maxTime,'.pdf',sep=''),dpi=600,width=10,height=8)
  
  #normalize each channel to -1, 1
  ndat = normalizeCh(origdat)
  
  #plot individual trends
  out = indivFitAndPlot(ndat,featShort,maxTime,logflag)
  indivFits = rbind(indivFits,out)
  
  #plot group trends
  dogIdx= ndat$subjid=="I004_A0001_D00" | ndat$subjid=="I004_A0002_D00" | ndat$subjid=="I004_A0003_D00" | ndat$subjid=="I004_A0004_D00"
  tmpdat = ndat[!dogIdx,]
  origdat = origdat[!dogIdx,]
  tmpdat_ch = summarySE(tmpdat, measurevar="feat", groupvars=c("time","subjid"))
  
  #plot mean across subjects after normalization
  tmpdat_subjch = summarySE(tmpdat_ch, measurevar="feat", groupvars=c("time"))
  feature = 'Line Length'
  xlabel='Time (Day)'
  ylabel=paste('Normalized ',featLabel,sep='')
  title=paste('Mean ',featLabel,' across channels, averaged across subjects',sep='')
  #g2 = ggplot(dat=tmpdat_subjch, aes(x=time, y=feat)) +
  #  geom_errorbar(aes(ymin=feat-se, ymax=feat+se), width=.1) +
  #  theme(legend.title=element_blank(),text = element_text(size=textsize),axis.text.x = element_text(angle=90, vjust=1)) + xlab(xlabel) + ylab(ylabel)  + scale_x_continuous(labels = formatter2) + ylim(0,1.5)
  #g2
  #ggsave(paste('MeanFeat-subj-',featShort,maxTime,'.pdf',sep=''),dpi=600,width=10,height=8)
  
  #plot mean across subjects with each subj
  xlabel='Time (Day)'
  ylabel=paste('Normalized ',featLabel,sep='')
  title=paste('Mean ',featLabel,' across channels',sep='')
  #g3 = ggplot(tmpdat_ch,aes(x=time,y=feat)) + geom_point(alpha=0.4) + geom_smooth(aes(colour=subjid),span=0.1,se=F) + scale_colour_discrete(guide=FALSE) + 
  #  theme(legend.title=element_blank(),text = element_text(size=textsize),axis.text.x = element_text(angle=90, vjust=1)) + xlab(xlabel) + ylab(ylabel)  + scale_x_continuous(labels = formatter2) + ylim(0,1.5)
  #g3
  #ggsave(paste('MeanFeat-subj-indiv-',featShort,maxTime,'.pdf',sep=''),dpi=600,width=10,height=8)
  tmpdat$feat_orig = origdat$feat
  return(tmpdat)
}

maxTime = 200
logflag = 0
idx = c(2,6)
fFeat = vector("list",length(featShortList))
allDat = data.frame(ch=integer(),
                    subjid=character())
for (i in 1:length(featShortList)){
  featShort = featShortList[i]
  print(paste('Processing',i,featShort),sep=' ')
  tmp = indivAndGroupPlot(featShort,featLabelList[i],maxTime,logflag)
  colnames(tmp)[idx] = c(featShort,paste(featShort,'_orig',sep=''))
  allDat = merge(allDat,tmp,all=T);
}

# ######## PLOT FEATURES###########

##### CORRELATIONS BETWEEN ALL FEATURES #####
cor(allDat[sapply(allDat, is.numeric)],method="pearson",use="complete.obs")


## PLOT COFVAR before normalization
allDatNorm = allDat[,c(1,2,3,4,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37)]
allDatOrig = allDat[,c(1,2,3,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38)]
#wide to long
allDatLong = reshape(allDatNorm, varying = featShortList, v.names = "value",
                     timevar = "feature", times = featShortList, direction = "long")

allDatLongOrig = reshape(allDatOrig, varying = paste(featShortList,'_orig',sep=''), v.names = "value",
                         timevar = "feature", times = paste(featShortList,'_orig',sep=''), direction = "long")


###### Cofvar #####
MAXTIMEPLOT = 500
tmp = ddply(allDatLongOrig,.(subjid,time,feature),summarize,meanfeat = mean(value),sd = sd(value), cofvar = sd/meanfeat)
origfeatlist = featLabelList;
names(origfeatlist) = paste(featShortList,'_orig',sep='')
tmp$feature = revalue(tmp$feature,replace=origfeatlist)
tmp$feature = factor(tmp$feature,levels=featLabelList)

#### normalize cofvar
subjects = unique(tmp$subjid)
features = unique(tmp$feature)
for (i in 1:length(subjects)) {
  subj = subjects[i]
  feat = tmp[tmp$subjid==subj,]
  for (f in 1:length(features)) {
    tmpfeat = feat[feat$feature == features[f],'cofvar']
    q = quantile(tmpfeat,c(0.025,0.975),na.rm=TRUE)
    tmp[tmp$subjid==subj & tmp$feature==features[f],'cofvar'] = (tmpfeat-q[1])/(q[2]-q[1])
  }
}

#get threshold
tmpsubj = ddply(tmp,.(time,feature),summarize,meanfeat = median(cofvar,na.rm=T))
uniqueFeat = unique(tmpsubj$feature)
tmp$upper_thres = NA
meannullD = NA
sdnullD = NA
for (i in 1:length(uniqueFeat)) {
  testdat = tmpsubj[tmpsubj$feature==uniqueFeat[i],]
  nullD = testdat[testdat$time>200 & testdat$time<MAXTIMEPLOT,3]
  meannullD[i] = mean(nullD,na.rm=T)
  sdnullD[i] = sd(nullD,na.rm=T)
  upper_thres = meannullD[i] + 1.645*sdnullD[i]
  lower_thres =  meannullD[i] - 1.645*sdnullD[i]
  tmp[tmp$feature==uniqueFeat[i],'upper_thres'] = upper_thres
  tmp[tmp$feature==uniqueFeat[i],'lower_thres'] = lower_thres
}
covnullD = data.frame(feature=uniqueFeat,mean=meannullD,sd=sdnullD)

mean_se <- function(x, mult = 1) {  
  x <- na.omit(x)
  se <- mult * sqrt(var(x) / length(x))
  mean <- mean(x)
  data.frame(y = mean, ymin = mean - se, ymax = mean + se)
}


g = ggplot(tmp[(tmp$feature %in% timeFeats),],aes(x=time,y=cofvar))  + stat_smooth(aes(colour="Individual Trends",group=subjid),method="loess",span=smoothSpan,se=F,size=0.4,alpha=0.6,level=0.95)    + stat_summary(fun.data=mean_se,geom="ribbon",alpha=0.35)   + stat_smooth(aes(colour="Mean"),se=F,alpha=0.8,level=0.95)  + facet_wrap(~feature) +
  geom_hline(aes(colour="Reference Threshold",yintercept=upper_thres),size=1,linetype="dotted")  + theme(plot.title=element_text(size=18),axis.text = element_text(size=12),axis.title = element_text(size=14),strip.text.x = element_text(size = 14),legend.position='top',legend.text=element_text(size=12)) + coord_cartesian(ylim = c(0,1.25),xlim=c(0,MAXTIMEPLOT)) +
  xlab('Time (Days)') + ylab('Normalized Coefficient of Variation') + ggtitle('Spatial Variation - Time Features') + scale_colour_manual(name="",values=c("Individual Trends"="orange","Mean" = "blue","Reference Threshold"="red"))
g = g + guides(colour = guide_legend(override.aes = list(linetype = c(1,1,3))))
ggsave(paste('alltime_cofvar',MAXTIMEPLOT,'.pdf',sep=''),dpi=600,height=6,width=8)



timeToStability = NULL
features = timeFeats
fits = ggplot_build(g)$data[[3]] #get smoothed fits
for (i in 1:length(features)) {
  thres = tmp[tmp$feature==features[i],'upper_thres'][1] #get upper thres
  panelfits = fits[fits$PANEL==i,]
  ## find flattening (diff = 0)
  l0 = which(diff(panelfits$y)<0)
  g0 = which(diff(panelfits$y)>0)
  inflexion = min(l0[which(diff(l0)>1)[1]],g0[which(diff(g0)>1)[1]],na.rm=TRUE)
  dtix= panelfits$x[inflexion]
  dtiy = panelfits$y[inflexion]
  tmpfits = panelfits[panelfits$y>thres & panelfits$x<MAXTIMEPLOT,] #get first fits
  xy = panelfits[c(dim(tmpfits)[1],dim(tmpfits)[1]+1),2:3] #calc slope
  slope = (xy[1,2] - xy[2,2])/(xy[2,1]-xy[1,1])
  cross = xy[1,1] + (xy[1,2] - thres) / slope #interpolate
  tmptts= data.frame(feature=features[i],tts = cross,t1 = tmpfits[1,]$y,td1 = panelfits$y[1],tc = thres,dtix=dtix,dtiy=dtiy, type='Cofvar')
  timeToStability = rbind(timeToStability,tmptts)
}

g = ggplot(tmp[(tmp$feature %in% derFeats[c(1,2)]),],aes(x=time,y=cofvar))   + stat_summary(fun.data=mean_se,geom="ribbon",alpha=0.35)     + facet_wrap(~feature) +
  geom_hline(aes(colour="Reference Threshold",yintercept=upper_thres),size=1,linetype="dotted")  + theme(plot.title=element_text(size=18),axis.text = element_text(size=12),axis.title = element_text(size=14),strip.text.x = element_text(size = 14),legend.position='top',legend.text=element_text(size=12)) + coord_cartesian(ylim = c(0,1.25),xlim=c(0,MAXTIMEPLOT)) +
  xlab('Time (Days)') + ylab('Normalized Coefficient of Variation') + ggtitle('Spatial Variation - Misc Features') + scale_colour_manual(name="",values=c("Mean" = "blue","Reference Threshold"="red"))
g = g + guides(colour = guide_legend(override.aes = list(linetype = c(3))))
ggsave(paste('supp_derfeatscofvar',MAXTIMEPLOT,'.pdf',sep=''),dpi=600,height=4,width=8)

g = ggplot(tmp[(tmp$feature %in% derFeats[3]),],aes(x=time,y=cofvar))   + stat_smooth(aes(colour="Individual Trends",group=subjid),method="loess",span=smoothSpan,se=F,size=0.4,alpha=0.6,level=0.95)  + stat_summary(fun.data=mean_se,geom="ribbon",alpha=0.35)     + stat_smooth(aes(colour="Mean"),se=F,alpha=0.8,level=0.95)  +
  geom_hline(aes(colour="Reference Threshold",yintercept=upper_thres),size=1,linetype="dotted")  + theme(axis.text = element_text(size=12),axis.title = element_text(size=14),strip.text.x = element_text(size = 14),legend.position='top',legend.text=element_text(size=12)) + coord_cartesian(ylim = c(0,1.25),xlim=c(0,MAXTIMEPLOT)) +
  xlab('Time (Days)') + ylab('Normalized Coefficient of Variation') + ggtitle('Spatial Variation') + scale_colour_manual(name="",values=c("Individual Trends"="orange","Mean" = "blue","Reference Threshold"="red"))
g = g + guides(colour = guide_legend(override.aes = list(linetype = c(1,1,3))))

ggsave(paste('supp_rmsfeatscofvar',MAXTIMEPLOT,'.pdf',sep=''),dpi=600,height=4,width=4)


g = ggplot(tmp[(tmp$feature %in% derFeats),],aes(x=time,y=cofvar)) + stat_smooth(aes(colour="Individual Trends",group=subjid),method="loess",span=smoothSpan,se=F,size=0.4,alpha=0.6,level=0.95)    + stat_summary(fun.data=mean_se,geom="ribbon",alpha=0.35)     + stat_smooth(aes(colour="Mean"),se=F,alpha=0.8,level=0.95)  + facet_wrap(~feature) +
  geom_hline(aes(colour="Reference Threshold",yintercept=upper_thres),size=1,linetype="dotted")  + theme(axis.text = element_text(size=12),axis.title = element_text(size=14),strip.text.x = element_text(size = 14),legend.position='top',legend.text=element_text(size=12)) + coord_cartesian(ylim = c(0,1.25),xlim=c(0,MAXTIMEPLOT)) +
  xlab('Time (Days)') + ylab('Normalized Coefficient of Variation') + ggtitle('Spatial Variation') + scale_colour_manual(name="",values=c("Individual Trends"="orange","Mean" = "blue","Reference Threshold"="red"))
g 

features = derFeats
fits = ggplot_build(g)$data[[3]] #get smoothed fits
for (i in 1:length(features)) {
  thres = tmp[tmp$feature==features[i],'upper_thres'][1] #get upper thres
  panelfits = fits[fits$PANEL==i,]
  ## find flattening (diff = 0)
  l0 = which(diff(panelfits$y)<0)
  g0 = which(diff(panelfits$y)>0)
  inflexion = min(l0[which(diff(l0)>1)[1]],g0[which(diff(g0)>1)[1]],na.rm=TRUE)
  dtix= panelfits$x[inflexion]
  dtiy = panelfits$y[inflexion]
  tmpfits = panelfits[panelfits$y>thres & panelfits$x<MAXTIMEPLOT,] #get first fits
  xy = panelfits[c(dim(tmpfits)[1],dim(tmpfits)[1]+1),2:3] #calc slope
  slope = (xy[1,2] - xy[2,2])/(xy[2,1]-xy[1,1])
  cross = xy[1,1] + (xy[1,2] - thres) / slope #interpolate
  tmptts= data.frame(feature=features[i],tts = cross,t1 = tmpfits[1,]$y,td1 = panelfits$y[1],tc = thres,dtix=dtix,dtiy=dtiy, type='Cofvar')
  timeToStability = rbind(timeToStability,tmptts)
}



g = ggplot(tmp[(tmp$feature %in% powerFeats),],aes(x=time,y=cofvar)) + stat_smooth(aes(colour="Individual Trends",group=subjid),method="loess",span=smoothSpan,se=F,size=0.4,alpha=0.6,level=0.95)    + stat_summary(fun.data=mean_se,geom="ribbon",alpha=0.35)    + stat_smooth(aes(colour="Mean"),se=F,alpha=0.8,level=0.95)  + facet_wrap(~feature) +
  geom_hline(aes(colour="Reference Threshold",yintercept=upper_thres),size=1,linetype="dotted")  + theme(legend.title=element_text(size=14),axis.text = element_text(size=12),axis.title = element_text(size=14),strip.text.x = element_text(size = 14),legend.position='top',legend.text=element_text(size=12)) + coord_cartesian(ylim = c(0,1.25),xlim=c(0,MAXTIMEPLOT)) +
  xlab('Time (Days)') + ylab('Normalized Coefficient of Variation') + ggtitle('Spatial Variation - Spectral Features') + scale_colour_manual(name="",values=c("Individual Trends"="orange","Mean" = "blue","Reference Threshold"="red"))
g = g + guides(colour = guide_legend(override.aes = list(linetype = c(1,1,3))))
ggsave(paste('allpower_cofvar',MAXTIMEPLOT,'.pdf',sep=''),dpi=600,height=6,width=8)


features = powerFeats
fits = ggplot_build(g)$data[[3]] #get smoothed fits
for (i in 1:length(features)) {
  thres = tmp[tmp$feature==features[i],'upper_thres'][1] #get upper thres
  panelfits = fits[fits$PANEL==i,]
  ## find flattening (diff = 0)
  l0 = which(diff(panelfits$y)<0)
  g0 = which(diff(panelfits$y)>0)
  inflexion = min(l0[which(diff(l0)>1)[1]],g0[which(diff(g0)>1)[1]])
  dtix= panelfits$x[inflexion]
  dtiy = panelfits$y[inflexion]
  tmpfits = panelfits[panelfits$y>thres & panelfits$x<MAXTIMEPLOT,] #get first fits
  xy = panelfits[c(dim(tmpfits)[1],dim(tmpfits)[1]+1),2:3] #calc slope
  slope = (xy[1,2] - xy[2,2])/(xy[2,1]-xy[1,1])
  cross = xy[1,1] + (xy[1,2] - thres) / slope #interpolate
  tmptts= data.frame(feature=features[i],tts = cross,t1 = tmpfits[1,]$y,td1=panelfits$y[1],tc = thres,dtix=dtix,dtiy=dtiy,type='Cofvar')
  timeToStability = rbind(timeToStability,tmptts)
}

g = ggplot(tmp[(tmp$feature %in% power2Feats),],aes(x=time,y=cofvar)) + stat_smooth(aes(colour="Individual Trends",group=subjid),method="loess",span=smoothSpan,se=F,size=0.4,alpha=0.6,level=0.95)    + stat_summary(fun.data=mean_se,geom="ribbon",alpha=0.35)    + stat_smooth(aes(colour="Mean"),se=F,alpha=0.8,level=0.95)  + facet_wrap(~feature) +
  geom_hline(aes(colour="Reference Threshold",yintercept=upper_thres),size=1,linetype="dotted")  + theme(legend.title=element_text(size=14),axis.text = element_text(size=12),axis.title = element_text(size=14),strip.text.x = element_text(size = 14),legend.position='top',legend.text=element_text(size=12)) + coord_cartesian(ylim = c(0,1.25),xlim=c(0,MAXTIMEPLOT)) +
  xlab('Time (Days)') + ylab('Normalized Coefficient of Variation') + ggtitle('Spatial Variation') + scale_colour_manual(name="",values=c("Individual Trends"="orange","Mean" = "blue","Reference Threshold"="red"))
g
ggsave(paste('supp_allpower2_cofvar',MAXTIMEPLOT,'.pdf',sep=''),dpi=600,height=6,width=8)

features = power2Feats
fits = ggplot_build(g)$data[[3]] #get smoothed fits
for (i in 1:length(features)) {
  thres = tmp[tmp$feature==features[i],'upper_thres'][1] #get upper thres
  panelfits = fits[fits$PANEL==i,]
  ## find flattening (diff = 0)
  l0 = which(diff(panelfits$y)<0)
  g0 = which(diff(panelfits$y)>0)
  inflexion = min(l0[which(diff(l0)>1)[1]],g0[which(diff(g0)>1)[1]])
  dtix= panelfits$x[inflexion]
  dtiy = panelfits$y[inflexion]
  tmpfits = panelfits[panelfits$y>thres & panelfits$x<MAXTIMEPLOT,] #get first fits
  xy = panelfits[c(dim(tmpfits)[1],dim(tmpfits)[1]+1),2:3] #calc slope
  slope = (xy[1,2] - xy[2,2])/(xy[2,1]-xy[1,1])
  cross = xy[1,1] + (xy[1,2] - thres) / slope #interpolate
  tmptts= data.frame(feature=features[i],tts = cross,t1 = tmpfits[1,]$y,td1=panelfits$y[1],tc = thres,dtix=dtix,dtiy=dtiy,type='Cofvar')
  timeToStability = rbind(timeToStability,tmptts)
}
#### Mean Feat ####
tmp = ddply(allDatLong,.(subjid,time,feature),summarize,meanfeat = mean(value),sd = sd(value))
featlist = featLabelList;
names(featlist) = featShortList
tmp$feature = revalue(tmp$feature,replace=featlist)
tmp$feature = factor(tmp$feature,levels=featLabelList)

#get threshold
tmpsubj = ddply(tmp,.(time,feature),summarize,meanfeat = median(meanfeat,na.rm=T))
uniqueFeat = unique(tmpsubj$feature)
tmp$upper_thres = NA
meannullD = NA
sdnullD = NA
for (i in 1:length(uniqueFeat)) {
  testdat = tmpsubj[tmpsubj$feature==uniqueFeat[i],]
  nullD = testdat[testdat$time>200 & testdat$time<MAXTIMEPLOT,3]
  meannullD[i] = mean(nullD,na.rm=T)
  sdnullD[i] = sd(nullD,na.rm=T)
  upper_thres = 1.645*sdnullD[i] + meannullD[i]
  tmp[tmp$feature==uniqueFeat[i],'upper_thres'] = upper_thres
}
featnullD = data.frame(feature=uniqueFeat,mean=meannullD,sd=sdnullD)


g = ggplot(tmp[(tmp$feature %in% timeFeats),],aes(x=time,y=meanfeat)) + stat_smooth(aes(colour="Individual Trends",group=subjid),method="loess",span=smoothSpan,se=F,size=0.4,alpha=0.6,level=0.95)   + stat_summary(fun.data=mean_se,geom="ribbon",alpha=0.35)    + stat_smooth(aes(colour="Mean"),se=F,alpha=0.8,level=0.95)  + facet_wrap(~feature) +
  geom_hline(aes(colour="Reference Threshold",yintercept=upper_thres),size=1,linetype="dotted")  + theme(axis.text = element_text(size=12),axis.title = element_text(size=14),strip.text.x = element_text(size = 14),legend.position='top',legend.text=element_text(size=12))  + coord_cartesian(ylim = c(0.15,1),xlim=c(0,MAXTIMEPLOT)) +
  xlab('Time (Days)') + ylab('Normalized Feature Value') + ggtitle('Mean value - Time Features') + scale_colour_manual(name="",values=c("Individual Trends"="orange","Mean" = "blue","Reference Threshold"="red"))
g = g + guides(colour = guide_legend(override.aes = list(linetype = c(1,1,3))))

ggsave(paste('alltime_meanfeat',MAXTIMEPLOT,'.pdf',sep=''),dpi=600,height=6,width=8)



features = timeFeats
fits = ggplot_build(g)$data[[3]] #get smoothed fits
for (i in 1:length(features)) {
  thres = tmp[tmp$feature==features[i],'upper_thres'][1] #get upper thres
  panelfits = fits[fits$PANEL==i,]
  ## find flattening (diff = 0)
  l0 = which(diff(panelfits$y)<0)
  g0 = which(diff(panelfits$y)>0)
  inflexion = min(l0[which(diff(l0)>1)[1]],g0[which(diff(g0)>1)[1]])
  dtix= panelfits$x[inflexion]
  dtiy = panelfits$y[inflexion]
  tmpfits = panelfits[panelfits$y>thres & panelfits$x<MAXTIMEPLOT,] #get first fits
  xy = panelfits[c(dim(tmpfits)[1],dim(tmpfits)[1]+1),2:3] #calc slope
  slope = (xy[1,2] - xy[2,2])/(xy[2,1]-xy[1,1])
  cross = xy[1,1] + (xy[1,2] - thres) / slope #interpolate
  tmptts= data.frame(feature=features[i],tts = cross,t1 = tmpfits[1,]$y,td1=panelfits$y[1],tc = thres,dtix=dtix,dtiy=dtiy,type='feat')
  timeToStability = rbind(timeToStability,tmptts)
}

g = ggplot(tmp[(tmp$feature %in% derFeats[1:2]),],aes(x=time,y=meanfeat))  + stat_summary(fun.data=mean_se,geom="ribbon",alpha=0.35)     + facet_wrap(~feature) +
  geom_hline(aes(colour="Reference Threshold",yintercept=upper_thres),size=1,linetype="dotted")  + theme(axis.text = element_text(size=12),axis.title = element_text(size=14),strip.text.x = element_text(size = 14),legend.position='top',legend.text=element_text(size=12)) + coord_cartesian(ylim=c(0.4,0.6),xlim=c(0,MAXTIMEPLOT)) +
  xlab('Time (Days)') + ylab('Normalized Feature Value') + ggtitle('Mean value - Time Features') + scale_colour_manual(name="",values=c("Mean" = "blue","Reference Threshold"="red"))
g = g + guides(colour = guide_legend(override.aes = list(linetype = c(3))))

ggsave(paste('supp_derfeatsmean',MAXTIMEPLOT,'.pdf',sep=''),dpi=600,height=4,width=8)


g = ggplot(tmp[(tmp$feature %in% derFeats[3]),],aes(x=time,y=meanfeat)) + stat_smooth(aes(colour="Individual Trends",group=subjid),method="loess",span=smoothSpan,se=F,size=0.4,alpha=0.6,level=0.95)    + stat_summary(fun.data=mean_se,geom="ribbon",alpha=0.35)    + stat_smooth(aes(colour="Mean"),se=F,alpha=0.8,level=0.95)  +
  geom_hline(aes(colour="Reference Threshold",yintercept=upper_thres),size=1,linetype="dotted")  + theme(axis.text = element_text(size=12),axis.title = element_text(size=14),strip.text.x = element_text(size = 14),legend.position='top',legend.text=element_text(size=12)) + coord_cartesian(xlim=c(0,MAXTIMEPLOT)) +
  xlab('Time (Days)') + ylab('Normalized Feature Value') + ggtitle('Mean value - Time Features') + scale_colour_manual(name="",values=c("Individual Trends"="orange","Mean" = "blue","Reference Threshold"="red"))
g = g + guides(colour = guide_legend(override.aes = list(linetype = c(1,1,3))))

ggsave(paste('supp_rmsfeatsmean',MAXTIMEPLOT,'.pdf',sep=''),dpi=600,height=4,width=4)

g = ggplot(tmp[(tmp$feature %in% derFeats),],aes(x=time,y=meanfeat)) + stat_smooth(aes(colour="Individual Trends",group=subjid),method="loess",span=smoothSpan,se=F,size=0.4,alpha=0.6,level=0.95)    + stat_summary(fun.data=mean_se,geom="ribbon",alpha=0.35)    + stat_smooth(aes(colour="Mean"),se=F,alpha=0.8,level=0.95)  + facet_wrap(~feature) + 
  geom_hline(aes(colour="Reference Threshold",yintercept=upper_thres),size=1,linetype="dotted")  + theme(axis.text = element_text(size=12),axis.title = element_text(size=14),strip.text.x = element_text(size = 14),legend.position='top',legend.text=element_text(size=12)) + coord_cartesian(xlim=c(0,MAXTIMEPLOT)) +
  xlab('Time (Days)') + ylab('Normalized Feature Value') + ggtitle('Mean value - Time Features') + scale_colour_manual(name="",values=c("Individual Trends"="orange","Mean" = "blue","Reference Threshold"="red"))
g 

features = derFeats
fits = ggplot_build(g)$data[[3]] #get smoothed fits
for (i in 1:length(features)) {
  thres = tmp[tmp$feature==features[i],'upper_thres'][1] #get upper thres
  panelfits = fits[fits$PANEL==i,]
  ## find flattening (diff = 0)
  l0 = which(diff(panelfits$y)<0)
  g0 = which(diff(panelfits$y)>0)
  inflexion = min(l0[which(diff(l0)>1)[1]],g0[which(diff(g0)>1)[1]])
  dtix= panelfits$x[inflexion]
  dtiy = panelfits$y[inflexion]
  tmpfits = panelfits[panelfits$y>thres & panelfits$x<MAXTIMEPLOT,] #get first fits
  xy = panelfits[c(dim(tmpfits)[1],dim(tmpfits)[1]+1),2:3] #calc slope
  slope = (xy[1,2] - xy[2,2])/(xy[2,1]-xy[1,1])
  cross = xy[1,1] + (xy[1,2] - thres) / slope #interpolate
  tmptts= data.frame(feature=features[i],tts = cross,t1 = tmpfits[1,]$y,td1=panelfits$y[1],tc = thres,dtix=dtix,dtiy=dtiy,type='feat')
  timeToStability = rbind(timeToStability,tmptts)
}


g = ggplot(tmp[(tmp$feature %in% powerFeats),],aes(x=time,y=meanfeat)) + stat_smooth(aes(colour="Individual Trends",group=subjid),method="loess",span=smoothSpan,se=F,size=0.4,alpha=0.6,level=0.95)    + stat_summary(fun.data=mean_se,geom="ribbon",alpha=0.35)    + stat_smooth(aes(colour="Mean"),se=F,alpha=0.8,level=0.95)  + facet_wrap(~feature) +
  geom_hline(aes(colour="Reference Threshold",yintercept=upper_thres),size=1,linetype="dotted")  + theme(axis.text = element_text(size=12),axis.title = element_text(size=14),strip.text.x = element_text(size = 14),legend.position='top',legend.text=element_text(size=12)) + coord_cartesian(ylim = c(0.2,1),xlim=c(0,MAXTIMEPLOT)) +
  xlab('Time (Days)') + ylab('Normalized Feature Value') + ggtitle('Mean value - Spectral Features') + scale_colour_manual(name="",values=c("Individual Trends"="orange","Mean" = "blue","Reference Threshold"="red"))
g = g + guides(colour = guide_legend(override.aes = list(linetype = c(1,1,3))))
ggsave(paste('allpower_meanfeat',MAXTIMEPLOT,'.pdf',sep=''),dpi=600,height=6,width=8)




features = powerFeats
fits = ggplot_build(g)$data[[3]] #get smoothed fits
for (i in 1:length(features)) {
  thres = tmp[tmp$feature==features[i],'upper_thres'][1] #get upper thres
  panelfits = fits[fits$PANEL==i,]
  ## find flattening (diff = 0)
  l0 = which(diff(panelfits$y)<0)
  g0 = which(diff(panelfits$y)>0)
  inflexion = min(l0[which(diff(l0)>1)[1]],g0[which(diff(g0)>1)[1]])
  dtix= panelfits$x[inflexion]
  dtiy = panelfits$y[inflexion]
  tmpfits = panelfits[panelfits$y>thres & panelfits$x<MAXTIMEPLOT,] #get first fits
  xy = panelfits[c(dim(tmpfits)[1],dim(tmpfits)[1]+1),2:3] #calc slope
  slope = (xy[1,2] - xy[2,2])/(xy[2,1]-xy[1,1])
  cross = xy[1,1] + (xy[1,2] - thres) / slope #interpolate
  tmptts= data.frame(feature=features[i],tts = cross,t1 = tmpfits[1,]$y,td1=panelfits$y[1],tc = thres,dtix=dtix,dtiy=dtiy,type='feat')
  timeToStability = rbind(timeToStability,tmptts)
}

g = ggplot(tmp[(tmp$feature %in% power2Feats & tmp$time<500),],aes(x=time,y=meanfeat)) + stat_smooth(aes(colour="Individual Trends",group=subjid),method="loess",span=smoothSpan,se=F,size=0.4,alpha=0.6,level=0.95)    + stat_summary(fun.data=mean_se,geom="ribbon",alpha=0.35)    + stat_smooth(aes(colour="Mean"),se=F,alpha=0.8,level=0.95)  + facet_wrap(~feature) +
  geom_hline(aes(colour="Reference Threshold",yintercept=upper_thres),size=1,linetype="dotted")  + theme(axis.text = element_text(size=12),axis.title = element_text(size=14),strip.text.x = element_text(size = 14),legend.position='top',legend.text=element_text(size=12)) + coord_cartesian(ylim = c(0.2,.75),xlim=c(0,MAXTIMEPLOT)) +
  xlab('Time (Days)') + ylab('Normalized Feature Value') + ggtitle('Mean value - Spectral Features') + scale_colour_manual(name="",values=c("Individual Trends"="orange","Mean" = "blue","Reference Threshold"="red"))
g = g + guides(colour = guide_legend(override.aes = list(linetype = c(1,1,3))))
ggsave(paste('supp_allpower2_meanfeat',MAXTIMEPLOT,'.pdf',sep=''),dpi=600,height=6,width=8)

features = power2Feats
fits = ggplot_build(g)$data[[3]] #get smoothed fits
for (i in 1:length(features)) {
  thres = tmp[tmp$feature==features[i],'upper_thres'][1] #get upper thres
  panelfits = fits[fits$PANEL==i,]
  ## find flattening (diff = 0)
  l0 = which(diff(panelfits$y)<0)
  g0 = which(diff(panelfits$y)>0)
  inflexion = min(l0[which(diff(l0)>1)[1]],g0[which(diff(g0)>1)[1]])
  dtix= panelfits$x[inflexion]
  dtiy = panelfits$y[inflexion]
  tmpfits = panelfits[panelfits$y>thres & panelfits$x<MAXTIMEPLOT,] #get first fits
  xy = panelfits[c(dim(tmpfits)[1],dim(tmpfits)[1]+1),2:3] #calc slope
  slope = (xy[1,2] - xy[2,2])/(xy[2,1]-xy[1,1])
  cross = xy[1,1] + (xy[1,2] - thres) / slope #interpolate
  tmptts= data.frame(feature=features[i],tts = cross,t1 = tmpfits[1,]$y,td1=panelfits$y[1],tc = thres,dtix=dtix,dtiy=dtiy,type='feat')
  timeToStability = rbind(timeToStability,tmptts)
}



#### individual feature plots grid ####
indivFit= function(tmpdat,featName,logflag) {
  print(featName)
  tmpavg = ddply(tmpdat, .(subjid,time),summarise,feat=mean(get(featName)))
  tmpavg = groupedData(feat~time|subjid,data=tmpavg)
  if (logflag ==0 ) {
    mod = lmList(feat~time, data=tmpavg,na.action=na.omit)
  }
  if (logflag ==1 ) {
    mod = lmList(log(feat)~time, data=tmpavg,na.action=na.omit)
  }
  if (logflag ==2 ) {
    mod = lmList(log(feat)~log(time), data=tmpavg,na.action=na.omit)
  }
  f1.lis = summary(mod)
  a = f1.lis$coefficients
  a = a[,,2]
  a = cbind(a,f1.lis$r.squared)
  colnames(a)[5] = 'r.squared'
  coeff = a[order(rownames(a)),]
  out = data.frame(subjid=rownames(coeff),feat=featName,b=coeff[,1],t=coeff[,3],p=coeff[,4],rsquared=coeff[,5])
  out$p.corrected = p.adjust(out$p,method='BH')
  return(out)
}
indivFitAndPlot = function(tmpdat,featName,maxTime,logflag) {
  tmpdat = tmpdat[tmpdat$time<maxTime,]
  indivFits = indivFit(tmpdat,featName,logflag)
  #plot fit
  tmpdatavg =  summarySE(tmpdat, measurevar=featName, groupvars=c("time","subjid"),na.rm=TRUE)
  #if (logflag ==0 ) {
  #  ggplot(tmpdatavg,aes(x=time,y=get(featName))) + geom_point() + geom_smooth(method='lm') + ylab(featName) + xlab('Time (Day)') + ggtitle(paste('Individual fits for ',featName,' by subject'))  + facet_wrap(~subjid,scales="free")
  #  ggsave(paste('indiv_fits_nolog',featName,'-',maxTime,'.pdf',sep=''),width=8,height=10,dpi=600)  
 # }
  #if (logflag ==1 ) {
  #  ggplot(tmpdatavg,aes(x=time,y=log(et(featName)))) + geom_point() + geom_smooth(method='lm') + ylab(paste('log(',featName,')')) + xlab('Time (Day)') + ggtitle(paste('Individual fits for ',featName,' by subject'))  + facet_wrap(~subjid,scales="free")
  #  ggsave(paste('indiv_fits_log',featName,'-',maxTime,'.pdf',sep=''),width=8,height=10,dpi=600)  
  #}
 # if (logflag ==2 ) {
  #  ggplot(tmpdatavg,aes(x=log(time),y=log(et(featName)))) + geom_point() + geom_smooth(method='lm') + ylab(paste('log(',featName,')')) + xlab('log(Time (Day))') + ggtitle(paste('Individual fits for ',featName,' by subject'))  + facet_wrap(~subjid,scales="free")
  #  ggsave(paste('indiv_fits_loglog',featName,'-',maxTime,'.pdf',sep=''),width=8,height=10,dpi=600)
  #}
  head(indivFits)
  return(indivFits)
}

summ = NULL
# load 
maxTime = 75
id = '24hr2'
logflag = 0
indivFitsAll = NULL

for (i in 1:length(featShortList)) {
  print(paste(i,featShortList[i]))
  featName = featShortList[i]
  out = indivFitAndPlot(allDat,featName,maxTime,logflag)
  head(out)
  indivFitsAll = rbind(indivFitsAll,out)
}

indivFitsAll$p.corrected = p.adjust(indivFitsAll$p,method='BH')
indivFitsAll$isSig = indivFitsAll$p.corrected<=0.05
indivFitsAll$'p<0.05' = 'black'
indivFitsAll$'p<0.05'[!indivFitsAll$isSig] = NA
indivFitsAll$onlySigBetas = indivFitsAll$b
indivFitsAll$onlySigBetas[!indivFitsAll$isSig] = NA
write.csv(indivFitsAll,paste('indivFitsAll-',id,logflag,maxTime,'.csv',sep=''))

indivFitsAll$subjid = revalue(indivFitsAll$subjid,c("NVC1001_23_002" = "1","NVC1001_23_003" = "2","NVC1001_23_004" = "3","NVC1001_23_005" = "4","NVC1001_23_006" = "5",
                                              "NVC1001_23_007" = "6","NVC1001_24_001" = "7","NVC1001_24_002" = "8","NVC1001_24_004" = "9","NVC1001_24_005" = "10",
                                              "NVC1001_25_001" = "11","NVC1001_25_002" = "12","NVC1001_25_003" = "13","NVC1001_25_004" = "14","NVC1001_25_005" = "15"))

# ggplot(data =  indivFitsAll, aes(x = feat, y = subjid,fill=b)) +  geom_tile(data=indivFitsAll,aes(colour='p<0.05')) +
#   geom_tile(dat=indivFitsAll,color=indivFitsAll$'p<0.05',size=0.8,height=0.95,width=0.94) + scale_colour_discrete(name='Outline') +
#   scale_fill_gradientn(name='b',colours=c("blue","cyan","white", "yellow","red"), 
#                        values=rescale(c(-1,0-.Machine$double.eps,0,0+.Machine$double.eps,1)),limits=c(-0.015,0.015)) + 
#   theme(legend.text=element_text(size=12),legend.title=element_text(size=12),axis.title.y = element_text(size=14), axis.text.x = element_text(angle=45, vjust=0.5,size=14), axis.text.y = element_text(size=14)) + ylab('Patient') + xlab('Feature') +
#   ggtitle('Rate of feature change in initial 75 days')

ggplot(data =  indivFitsAll, aes(x = feat, y = subjid,fill=onlySigBetas)) +  geom_tile(colour='white') +
  scale_fill_gradientn(name='b',colours=c("blue","cyan","white", "yellow","red"), 
                       values=rescale(c(-1,0-.Machine$double.eps,0,0+.Machine$double.eps,1)),limits=c(-0.015,0.015)) + 
  theme(plot.title = element_text(size = 16),legend.text=element_text(size=12),legend.title=element_text(size=12),axis.title.y = element_text(size=14), axis.title.x = element_text(size=14),axis.text.x = element_text(angle=45, vjust=0.5,size=14), axis.text.y = element_text(size=14)) + ylab('Patient') + xlab('Feature') +
  ggtitle('Rate of feature change in initial 75 days')

ggsave(paste('indivFitsAll-heatmap-log',id,logflag,maxTime,'.pdf',sep=''),width=8,height=11,dpi=600)


ggplot(data =  indivFitsAll, aes(x = feat, y = subjid,fill=onlySigBetas)) +  geom_tile(colour='white') + geom_text(aes(label = round(b, 3))) +
  scale_fill_gradientn(name='b',colours=c("blue","cyan","white", "yellow","red"), 
                       values=rescale(c(-1,0-.Machine$double.eps,0,0+.Machine$double.eps,1)),limits=c(-0.015,0.015)) + 
  theme(plot.title = element_text(size = 16),legend.text=element_text(size=12),legend.title=element_text(size=12),axis.title.y = element_text(size=14), axis.title.x = element_text(size=14),axis.text.x = element_text(angle=45, vjust=0.5,size=14), axis.text.y = element_text(size=14)) + ylab('Patient') + xlab('Feature') +
  ggtitle('Rate of feature change in initial 75 days with beta values')

ggsave(paste('indivFitsAll-heatmapb-log',id,logflag,maxTime,'.pdf',sep=''),width=8,height=11,dpi=600)


ggplot(data =  indivFitsAll, aes(x = feat, y = subjid,fill=onlySigBetas)) +  geom_tile(colour='white') + geom_text(aes(label = round(t, 3))) +
  scale_fill_gradientn(name='b',colours=c("blue","cyan","white", "yellow","red"), 
                       values=rescale(c(-1,0-.Machine$double.eps,0,0+.Machine$double.eps,1)),limits=c(-0.015,0.015)) + 
  theme(plot.title = element_text(size = 16),legend.text=element_text(size=12),legend.title=element_text(size=12),axis.title.y = element_text(size=14), axis.title.x = element_text(size=14),axis.text.x = element_text(angle=45, vjust=0.5,size=14), axis.text.y = element_text(size=14)) + ylab('Patient') + xlab('Feature') +
  ggtitle('Rate of feature change in initial 75 days with t statistic')

ggsave(paste('indivFitsAll-heatmapt-log',id,logflag,maxTime,'.pdf',sep=''),width=8,height=11,dpi=600)

indivFitsAll$p.corrected = round(indivFitsAll$p.correct,3)
indivFitsAll$p.corrected[indivFitsAll$p.corrected < 0.001] = '<0.001'
ggplot(data =  indivFitsAll, aes(x = feat, y = subjid,fill=onlySigBetas)) +  geom_tile(colour='white') + geom_text(aes(label = p.corrected)) +
  scale_fill_gradientn(name='b',colours=c("blue","cyan","white", "yellow","red"), 
                       values=rescale(c(-1,0-.Machine$double.eps,0,0+.Machine$double.eps,1)),limits=c(-0.015,0.015)) + 
  theme(plot.title = element_text(size = 16),legend.text=element_text(size=12),legend.title=element_text(size=12),axis.title.y = element_text(size=14), axis.title.x = element_text(size=14),axis.text.x = element_text(angle=45, vjust=0.5,size=14), axis.text.y = element_text(size=14)) + ylab('Patient') + xlab('Feature') +
  ggtitle('Rate of feature change in initial 75 days with corrected p values')

ggsave(paste('indivFitsAll-heatmapp-log',id,logflag,maxTime,'.pdf',sep=''),width=8,height=11,dpi=600)

#### INDIVIDUAL FEATURE PLOTS AND DAYS TO STABILIZATION #####
MAXTIMEPLOT=500
maxTime =500
smoothSpan = 0.65
calc_thres = function(tmp) {
  nullD = tmp[tmp$time>200,]
  if (nrow(nullD) > 50){
    #ggplot(data=nullD,aes(x=time,y=feat)) + geom_point() + geom_smooth(method=lm)
    mod = lm(data=nullD,feat~time)
    adj_r = summary(mod)$adj.r.squared
    modsumm = summary(mod)$coefficients
    # lin.mod = lm(data=tmp,feat~time)
    # smod = segmented(lin.mod,seg.Z = ~time, psi=14)
    # plot(data=tmp,feat~time)
    # plot(smod, add=T)
    tmpsumm = summary(mod)
    p_val = tmpsumm$coefficients[8]
    b_val = tmpsumm$coefficients[2]
    mod2 = lm(data=tmp,feat~time)
    mod2summ = summary(mod2)$coefficients
    b = c(modsumm[2] - modsumm[4],modsumm[2] + modsumm[4])
    b2 = c(mod2summ[2] - mod2summ[4],mod2summ[2] + mod2summ[4])
  } else {
    p_val = -999
    b_val = -999
    adj_r = -999
  }
  meannullD = mean(nullD[,3],na.rm=T)
  sdnullD = sd(nullD[,3],na.rm=T)
  upper_thres = meannullD + 1.645*sdnullD
  lower_thres =  meannullD - 1.645*sdnullD
  return(c(upper_thres,lower_thres,p_val,b_val,adj_r))
}
getStability= function(dat, featName,uniqueSubj) {

  tmpdat_ch = summarySE(dat, measurevar="feat", groupvars=c("time","subjid"))
  tmpdat_ch$subjid = revalue(tmpdat_ch$subjid,c("NVC1001_23_002" = "1","NVC1001_23_003" = "2","NVC1001_23_004" = "3","NVC1001_23_005" = "4","NVC1001_23_006" = "5",
                                                "NVC1001_23_007" = "6","NVC1001_24_001" = "7","NVC1001_24_002" = "8","NVC1001_24_004" = "9","NVC1001_24_005" = "10",
                                                "NVC1001_25_001" = "11","NVC1001_25_002" = "12","NVC1001_25_003" = "13","NVC1001_25_004" = "14","NVC1001_25_005" = "15"))
  
  #mean value per feat, per subject
  colnames(tmpdat_ch)[4]=featName
  tmpavg = ddply(tmpdat_ch, .(subjid,time),summarise,feat=mean(get(featName),na.rm=TRUE))
  tmpavg = groupedData(feat~time|subjid,data=tmpavg)
  gthres = gapply(tmpavg,FUN=calc_thres)
  gthres=data.frame(gthres)
  gthres$subjid = rownames(gthres)
  #g = ggplot(tmpavg,aes(x=time,y=feat)) + geom_point() + stat_smooth(aes(colour="Smoothed"),method='loess',se=F,alpha=0.8,level=0.95) + 
  #  geom_hline(dat=gthres,aes(colour="Reference Threshold",yintercept=gthres),size=1,linetype="dotted") + scale_colour_manual(name="",values=c("Smoothed" = "blue","Reference Threshold"="red")) + facet_wrap(~subjid)
  #g  
  #ggsave(paste(featName,'-daystostabilization',maxTime,'.pdf',sep=''),width=8,height=10,dpi=600) 

  timeToStability = NULL
  for (tsubj in 1:15)
  {
    tmp = tmpavg[tmpavg$subjid==tsubj,]
    tmp_thres = calc_thres(tmp)
    upper_thres = tmp_thres[1]
    lower_thres = tmp_thres[2]
    p_val = tmp_thres[3]
    b_val = tmp_thres[4]
    adj_r = tmp_thres[5]
    g = ggplot(tmp,aes(x=time,y=feat)) + geom_point() + stat_smooth(aes(colour="Mean"),method='loess',span=smoothSpan,se=F,alpha=0.8,level=0.95) +
      geom_hline(aes(colour="Reference Threshold",yintercept=upper_thres),size=1,linetype="dotted")  + 
      geom_hline(aes(colour="Reference Threshold",yintercept=lower_thres),size=1,linetype="dotted") +
      theme(plot.title=element_text(size=18),axis.text = element_text(size=12),axis.title = element_text(size=14),strip.text.x = element_text(size = 14),legend.position='top',legend.text=element_text(size=12)) + coord_cartesian(xlim=c(0,MAXTIMEPLOT)) +
      xlab('Time (Days)') + ylab('Normalized Value') + ggtitle(paste('Subj ',uniqueSubj[i], featName)) + scale_colour_manual(name="",values=c("Individual Trends"="orange","Mean" = "blue","Reference Threshold"="red"))
    g
    #ggsave(paste(uniqueSubj[i],featName,'-',maxTime,'-smoothed.pdf',sep=''),width=8,height=10,dpi=600) 
    fits = ggplot_build(g)$data[[2]] #get smoothed fits
    panelfits = fits[fits$PANEL==1,]
    ## find flattening (diff = 0)
    a = diff(panelfits$y)
    inflexion = min(which(a[1:length(a)-1]<0 & a[2:length(a)]>0), which(a[1:length(a)-1]>0 & a[2:length(a)]<0),na.rm=TRUE)
    dtix= panelfits$x[inflexion]
    dtiy = panelfits$y[inflexion]
    #g = g + geom_vline(aes(xintercept=dtix))
    #ggsave(paste(uniqueSubj[i],featName,'-',maxTime,'-indivinflexion.pdf',sep=''),width=8,height=10,dpi=600) 
    cross = NA
    if (!is.na(upper_thres)) {
      if (fits$y[1] > upper_thres){
        tmpfits = panelfits[panelfits$y>upper_thres & panelfits$x<MAXTIMEPLOT,] #get first fits
        xy = panelfits[c(dim(tmpfits)[1],dim(tmpfits)[1]+1),2:3] #calc slope
        slope = (xy[1,2] - xy[2,2])/(xy[2,1]-xy[1,1])
        cross = xy[1,1] + (xy[1,2] - upper_thres) / slope #interpolate
      }else{
        tmpfits = panelfits[panelfits$y<lower_thres & panelfits$x<MAXTIMEPLOT,] #get first fits
        xy = panelfits[c(dim(tmpfits)[1],dim(tmpfits)[1]+1),2:3] #calc slope
        slope = (xy[2,2] - xy[1,2])/(xy[2,1]-xy[1,1])
        cross = xy[1,1] + (xy[2,2] - lower_thres) / slope #interpolate
      }
    }
    if (tsubj %in% uniqueSubj) {
      #if (tsubj==6 && featName=='LL'){
       # dtix = 200
      #}
      #if (tsubj==10 && featName=='LL'){
      #  dtix = 220
      #}
      #if (tsubj==7 && featName=='Gamma'){
      #  dtix = 75
      #}
      if (p_val == -999 | (p_val < 0.05 & adj_r>0.5) )  {
        tmptts= data.frame(feature=featName,tts = NA,t1 = tmpfits[1,]$y,td1 = panelfits$y[1],tc = NA,p=p_val<0.05,beta=b_val,adjr=adj_r,tcl=NA,dtix=NA,dtiy=dtiy, type='feat')
      } else {
        tmptts= data.frame(feature=featName,tts = cross,t1 = tmpfits[1,]$y,td1 = panelfits$y[1],tc = upper_thres,p=p_val<0.05,beta=b_val,adjr=adj_r,tcl=lower_thres,dtix=dtix,dtiy=dtiy, type='feat')
      }
    }else {
      if (p_val == -999 | (p_val < 0.05 & adj_r>0.5) ) {#& abs(b_val)>1/300
        tmptts= data.frame(feature=featName,tts = NA,t1 = tmpfits[1,]$y,td1 = panelfits$y[1],tc = NA,p=p_val,beta=b_val,adjr=adj_r,tcl=NA,dtix=NA,dtiy=dtiy, type='feat')
      } else {
        tmptts= data.frame(feature=featName,tts = cross,t1 = tmpfits[1,]$y,td1 = panelfits$y[1],tc = upper_thres,p=p_val,beta=b_val,adjr=adj_r,tcl=lower_thres,dtix=NA,dtiy=dtiy, type='feat')
      }
    }
    timeToStability = rbind(timeToStability,tmptts)
  }
  return(timeToStability)
}

indivStabilityandPlot = function(featName,featLabel,stability,uniqueSubj) {
  id = '24hr2'
  dat = as.data.frame(readMat(paste("combined_",featName,"_",id,".mat",sep='')))
  #remove outliers
  dat = removeOutliers(dat)
  dat = normalizeCh(dat)
  ttstability = getStability(dat,featName,uniqueSubj)
  ttstability$subjid = c(1:15)
  stability = rbind(stability,ttstability)
  #plot fit
  tmpdatavg =  summarySE(dat, measurevar="feat", groupvars=c("time","subjid"),na.rm=TRUE)
  tmpdatavg$subjid = revalue(tmpdatavg$subjid,c("NVC1001_23_002" = "1","NVC1001_23_003" = "2","NVC1001_23_004" = "3","NVC1001_23_005" = "4","NVC1001_23_006" = "5",
                                                "NVC1001_23_007" = "6","NVC1001_24_001" = "7","NVC1001_24_002" = "8","NVC1001_24_004" = "9","NVC1001_24_005" = "10",
                                                "NVC1001_25_001" = "11","NVC1001_25_002" = "12","NVC1001_25_003" = "13","NVC1001_25_004" = "14","NVC1001_25_005" = "15"))
  merged = merge(tmpdatavg,ttstability,by="subjid",all=TRUE)
  g = ggplot(merged,aes(x=time,y=feat)) + geom_point() + stat_smooth(aes(colour="Smoothed Mean"),method="loess",span=smoothSpan,se=F,alpha=0.8,level=0.95)  + 
    geom_hline(aes(colour="Reference Threshold",yintercept=tc),size=1.25,linetype="dotted") +
    geom_hline(aes(colour="Reference Threshold",yintercept=tcl),size=1.25,linetype="dotted") + 
    geom_vline(aes(colour='Reference Point',xintercept=tts),size=1.25) +
    geom_vline(aes(colour='Inflexion Point',xintercept=dtix),size=1.25) + 
    facet_wrap(~subjid,scales="free_y") + coord_cartesian(xlim=c(0,MAXTIMEPLOT))  + 
    ylab(paste('Normalized ',featLabel, sep='')) + xlab('Time (Days)') + ggtitle(paste('Individual fits for ',featLabel,' by subject')) +
    scale_colour_manual(name="",values=c("Smoothed Mean" = "blue","Reference Threshold"="red","Reference Point"="orange","Inflexion Point"="#009E73"))
  ggsave(paste('allsubj_labeled',featLabel,'-',maxTime,'.pdf',sep=''),width=10,height=10,dpi=600)  
  return(stability)
}

featShortList2 = featShortList[c(1,2,3,4,5,6,7,8,9,12)]
featLabelList2 = featLabelList[c(1,2,3,4,5,6,7,8,9,12)]
stability = NULL
smoothSpanList = c(0.70,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65)
uniqueSubjList = list()
uniqueSubjList[[1]] = 1:15#c(5,6,7,8,9,10,11,13)
uniqueSubjList[[2]] = 1:15#c(7,9,14)
uniqueSubjList[[3]] = 1:15#c(4,6,7,8,9,11,13,14)
uniqueSubjList[[4]] = 1:15#c(2,6,7,8,9,13,14,15)
uniqueSubjList[[5]] = 1:15#c(4,6,7,8,9,11,13,14,15)
uniqueSubjList[[6]] = 1:15#c(2,4,6,7,8,9,11,13,14,15)
uniqueSubjList[[7]] = 1:15#c(2,4,6,7,8,9,11,13,14,15)
uniqueSubjList[[8]] = 1:15#c(2,6,7,8,9,10,11,12,13,14)
uniqueSubjList[[9]] = 1:15# c(3,6,7,8,9,11,13)
uniqueSubjList[[10]] = 1:15#c(1,4,7,8,9,11,13)

for (i in 1:length(featShortList2)){
  smoothSpan = smoothSpanList[i]
  uniqueSubj = uniqueSubjList[[i]]
  featName = featShortList2[i]
  featLabel = featLabelList2[i]
  print(paste(featName,featLabel))
  stability = indivStabilityandPlot(featName,featLabel,stability,uniqueSubj)
}

ggplot(data =  stability, aes(x = feature, y = subjid)) +  geom_tile(colour='white') + geom_text(aes(label = round(tts,3))) + ylab('Patient') + xlab('Feature') + ggtitle('Days to stability')
ggsave(paste('indiv-DaysToStability-feat',id,logflag,maxTime,'.pdf',sep=''),width=8,height=11,dpi=600)

ggplot(data =  stability, aes(x = feature, y = subjid)) +  geom_tile(colour='white') + geom_text(aes(label = round(dtix,3))) + ylab('Patient') + xlab('Feature') + ggtitle('Days to inflexion')
ggsave(paste('indiv-DaysToInflexion-feat',id,logflag,maxTime,'.pdf',sep=''),width=8,height=11,dpi=600)






####### SD FEAT ########
indivFit= function(tmpdat,featName,maxTime,logflag) {
  tmp = tmpdat[tmpdat$time<maxTime,]
  tmpavg = ddply(tmp, .(subjid,time),summarize,feat=mean(feat))
  tmpavg = groupedData(feat~time|subjid,data=tmpavg)
  if (logflag ==0 ) {
    mod = lmList(feat~time, data=tmpavg,na.action=na.omit)
  }
  if (logflag ==1 ) {
    mod = lmList(log(feat)~time, data=tmpavg,na.action=na.omit)
  }
  if (logflag ==2 ) {
    mod = lmList(log(feat)~log(time), data=tmpavg,na.action=na.omit)
  }
  f1.lis = summary(mod)
  a = f1.lis$coefficients
  a = a[,,2]
  a = cbind(a,f1.lis$r.squared)
  colnames(a)[5] = 'r.squared'
  coeff = a[order(rownames(a)),]
  out = data.frame(subjid=rownames(coeff),feat=featName,b=coeff[,1],t=coeff[,3],p=coeff[,4],rsquared=coeff[,5])
  return(out)
}
indivFitAndPlot = function(tmpdat,featName,maxTime,logflag) {
  indivFits = indivFit(tmpdat,featName,maxTime,logflag)
  #plot fit
  tmpdatavg =  summarySE(tmpdat, measurevar="feat", groupvars=c("time","subjid"),na.rm=TRUE)
  # if (logflag ==0 ) {
  #   ggplot(tmpdatavg,aes(x=time,y=feat)) + geom_point() + geom_smooth(method='lm') + ylab(featName) + xlab('Time (Day)') + ggtitle(paste('Individual fits for ',featName,' by subject')) + scale_x_continuous(labels = formatter2)  + facet_wrap(~subjid,scales="free")
  #   ggsave(paste('indiv_fits_nologsd',featName,'-',maxTime,'.pdf',sep=''),width=8,height=10,dpi=600)  
  # }
  # if (logflag ==1 ) {
  #   ggplot(tmpdatavg,aes(x=time,y=log(feat))) + geom_point() + geom_smooth(method='lm') + ylab(paste('log(',featName,')')) + xlab('Time (Day)') + ggtitle(paste('Individual fits for ',featName,' by subject')) + scale_x_continuous(labels = formatter2)  + facet_wrap(~subjid,scales="free")
  #   ggsave(paste('indiv_fits_logsd',featName,'-',maxTime,'.pdf',sep=''),width=8,height=10,dpi=600)  
  # }
  # if (logflag ==2 ) {
  #   ggplot(tmpdatavg,aes(x=log(time),y=log(feat))) + geom_point() + geom_smooth(method='lm') + ylab(paste('log(',featName,')')) + xlab('log(Time (Day))') + ggtitle(paste('Individual fits for ',featName,' by subject')) + scale_x_continuous(labels = formatter2)  + facet_wrap(~subjid,scales="free")
  #   ggsave(paste('indiv_fits_loglogsd',featName,'-',maxTime,'.pdf',sep=''),width=8,height=10,dpi=600)
  # }
  
  return(indivFits)
}
indivFits = NULL

id = '24hr2'
featShortList = c('ll','area','hwamp','energy','delta','theta','alpha','beta','gamma','higamma')
featLabelList = c('Line Length','Area','Halfwave Amp','Energy','Delta Power','Theta Power','Alpha Power','Beta Power','Gamma Power','High Gamma Power')
indivAndGroupPlotSD = function(featShort,featLabel,maxTime,logflag) {
  dat = as.data.frame(readMat(paste("combined_",featShort,"_",id,".mat",sep='')))
  #remove outliers
  dat[,2] = dat[,3]
  dat = dat[,-3]
  dat = removeOutliers(dat)
  origdat = dat
  dat = dat[dat$time<=maxTime,]
  textsize = 18
  #before normalization, calculate Coefficient of variation, since this needs to be on a ratio scale and nonnegative.
  tmpdat_ch = summarySE(dat, measurevar="feat", groupvars=c("time","subjid"))
  
  #plot Coeff of Var between channels, averaged across all patients 
  tmpdat_ch$cofvar = tmpdat_ch$sd/tmpdat_ch$feat 
  tmpdat_meancofvarch = summarySE(tmpdat_ch, measurevar="cofvar", groupvars=c("time"))
  xlabel='Time (Day)'
  ylabel=paste('Mean CoV - ',featLabel,sep='')
  title=paste('Mean CoV across channels, averaged across subjects',sep='')
  #g1 = ggplot(dat=tmpdat_meancofvarch, aes(x=time, y=cofvar)) + 
  #  geom_errorbar(aes(ymin=cofvar-se, ymax=cofvar+se), width=.1) +
  #  theme(legend.title=element_blank(),text = element_text(size=textsize),axis.text.x = element_text(angle=90, vjust=1)) + xlab(xlabel) + ylab(ylabel) + scale_x_continuous(labels = formatter2) 
  #g1
  #ggsave(paste('MeanCoV-subj-',featShort,maxTime,'.pdf',sep=''),dpi=600,width=10,height=8)
  
  #normalize each channel to -1, 1
  ndat = normalizeCh(dat)
  
  #plot individual trends
  out = indivFitAndPlot(ndat,featShort,maxTime,logflag)
  indivFits = rbind(indivFits,out)
  
  #plot group trends
  dogIdx= ndat$subjid=="I004_A0001_D00" | ndat$subjid=="I004_A0002_D00" | ndat$subjid=="I004_A0003_D00" | ndat$subjid=="I004_A0004_D00"
  tmpdat = ndat[!dogIdx,]
  origdat = origdat[!dogIdx,]
  tmpdat_ch = summarySE(tmpdat, measurevar="feat", groupvars=c("time","subjid"))
  
  #plot mean across subjects after normalization
  tmpdat_subjch = summarySE(tmpdat_ch, measurevar="feat", groupvars=c("time"))
  feature = 'Line Length'
  xlabel='Time (Day)'
  ylabel=paste('Normalized ',featLabel,sep='')
  title=paste('Mean ',featLabel,' across channels, averaged across subjects',sep='')
  #g2 = ggplot(dat=tmpdat_subjch, aes(x=time, y=feat)) +
  #  geom_errorbar(aes(ymin=feat-se, ymax=feat+se), width=.1) +
  #  theme(legend.title=element_blank(),text = element_text(size=textsize),axis.text.x = element_text(angle=90, vjust=1)) + xlab(xlabel) + ylab(ylabel)  + scale_x_continuous(labels = formatter2) + ylim(0,1.5)
  #g2
  #ggsave(paste('MeanFeat-subj-',featShort,maxTime,'.pdf',sep=''),dpi=600,width=10,height=8)
  
  #plot mean across subjects with each subj
  xlabel='Time (Day)'
  ylabel=paste('Normalized ',featLabel,sep='')
  title=paste('Mean ',featLabel,' across channels',sep='')
  #g3 = ggplot(tmpdat_ch,aes(x=time,y=feat)) + geom_point(alpha=0.4) + geom_smooth(aes(colour=subjid),span=0.1,se=F) + scale_colour_discrete(guide=FALSE) + 
  #  theme(legend.title=element_blank(),text = element_text(size=textsize),axis.text.x = element_text(angle=90, vjust=1)) + xlab(xlabel) + ylab(ylabel)  + scale_x_continuous(labels = formatter2) + ylim(0,1.5)
  #g3
  #ggsave(paste('MeanFeat-subj-indiv-',featShort,maxTime,'.pdf',sep=''),dpi=600,width=10,height=8)
  tmpdat$feat_orig = origdat$feat
  return(tmpdat)
}

maxTime = 2000
logflag = 2
fLL = indivAndGroupPlotSD(featShortList[1],featLabelList[1],maxTime,logflag)
fArea =indivAndGroupPlotSD(featShortList[2],featLabelList[2],maxTime,logflag)
fHW = indivAndGroupPlotSD(featShortList[3],featLabelList[3],maxTime,logflag)
fEnergy = indivAndGroupPlotSD(featShortList[4],featLabelList[4],maxTime,logflag)
fD = indivAndGroupPlotSD(featShortList[5],featLabelList[5],maxTime,logflag)
fT = indivAndGroupPlotSD(featShortList[6],featLabelList[6],maxTime,logflag)
fA = indivAndGroupPlotSD(featShortList[7],featLabelList[7],maxTime,logflag)
fB = indivAndGroupPlotSD(featShortList[8],featLabelList[8],maxTime,logflag)
fG = indivAndGroupPlotSD(featShortList[9],featLabelList[9],maxTime,logflag)
fHG = indivAndGroupPlotSD(featShortList[10],featLabelList[10],maxTime,logflag)


runGLM = function(dat) {
  featName = colnames(dat)[2]
  colnames(dat)[2] = 'feat'
  newdat = dat[dat$time<75,]
  dat = ddply(newdat, .(subjid,time),summarize,feat=mean(feat))
  dat = dat[!is.na(dat$feat),]
  tmpdat = dat
  #tmpdat$feat = log(dat$feat+1)
  # tmpdat$time = log(dat$time)
  humanAvgDatN = groupedData(feat~time|subjid,data=tmpdat ,labels=list(x="Time (Day)",y = "Feature"))
  mod = lme(feat ~ time,random=~time|subjid,data=humanAvgDatN)
  summ = summary(mod)$tTable
  rownames(summ) = c(paste('Intercept-',featName,sep=''),featName)
  return(summ)
}
#LMM
LLglm = runGLM(fLL)
Areaglm = runGLM(fArea)
HWglm = runGLM(fHW)
Energyglm = runGLM(fEnergy)
Dglm = runGLM(fD)
Tglm = runGLM(fT)
Aglm = runGLM(fA)
Bglm = runGLM(fB)
Gglm = runGLM(fG)
HGglm = runGLM(fHG)

tmp = rbind(LLglm,Areaglm,HWglm,Energyglm,Dglm,Tglm,Aglm,Bglm,Gglm,HGglm)
tmp = as.data.frame(tmp)
tmp$padjusted[seq(2,21,2)] = p.adjust(tmp[seq(2,21,2),'p-value'],'BH')
