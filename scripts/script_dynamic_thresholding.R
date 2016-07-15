#############################################
###############Whole DataSet#################
#############################################

require(arules)
require-(arulesViz)

load("~/Projects/MBA/data/MASTER_MATRIX.RData")
load("~/Projects/MBA/data/MASTER_MATRIX_PERMUTED.RData")

b.t<-as(as.data.frame(m.final),"transactions")
b.t.perm<-as(as.data.frame(m.final.perm),"transactions")
drug.idx<-grep("IC.50",labels(m.final)[[2]])
drugs<-labels(m.final)[[2]][drug.idx]
drugs.total<-c(paste(drugs,"Sensitive",sep="="),paste(drugs,"Resistant",sep="="))

rules<-apriori(b.t,parameter = list(support = (4/length(b.t)), confidence = 4/length(b.t), maxlen=2),appearance = list(rhs=drugs.total,default="lhs"))

temp.quality<-quality(rules)
temp.quality[,3]<-NA
temp.quality<-temp.quality[!duplicated(temp.quality), ]

for (i in 1:nrow(temp.quality)){
  print(i)
  rules.perm<-apriori(b.t.perm,parameter = list(support = temp.quality$support[i], confidence = temp.quality$confidence[i], maxlen=2))
  if (length(rules.perm)==0) next
  temp<-hist(log(quality(rules.perm)[,3]),breaks=80,plot=F)
  lift.estim.0.05<-exp(temp$breaks[which(cumsum(temp$counts)>sum(temp$counts)*0.95)[1]])
  temp.quality$lift[i]<-lift.estim.0.05
}
temp.quality$lift[temp.quality$lift<=1]<-NA
temp.quality$lift[is.na(temp.quality$lift)]<-min(temp.quality$lift,na.rm=T)


rm(b.t.perm)
rm(rules.perm)
gc()


status<-rep(NA,length(rules))
quality.rules<-quality(rules)
require(parallel)
require(arules)
cl<-makeCluster(4)
clusterExport(cl, "quality.rules")
clusterExport(cl, "temp.quality")
clusterExport(cl, "status")
for (i in 1:nrow(temp.quality)){
  print(i)
  clusterExport(cl, "i")
  temp.idx<-which(quality.rules$support==temp.quality$support[i] & quality.rules$confidence==temp.quality$confidence[i])
  clusterExport(cl, "temp.idx")
  status[temp.idx]<-parSapply (cl=cl, quality.rules$lift[temp.idx], function(lift) if(lift>temp.quality$lift[i]){temp.out<-"pass"}else{temp.out<-"fail"})
}
stopCluster(cl)
rm(cl)
gc()

rules.sig<-rules[which(status=="pass")]
save(rules.sig,file="~/Projects/MBA/data/RULES_SIGNIFICANT_Sup4_Conf4in689_DYNAMIC_THRESHOLD_FDR0.05.RData")
write(rules.sig,file="~/Projects/MBA/data/RULES_SIGNIFICANT_Sup4_Conf4in689_DYNAMIC_THRESHOLD_FDR0.05.tsv",sep='\t',row.names=F)
rm(rules)
gc()
#write.table(temp.quality,file="~/Projects/MBA/data/Thresholds_Sup_Conf_Lift_forFDR0.05_grid.tsv",sep='\t',row.names = F)








#############################################
#######Whole DataSet - Negative Control######
#############################################
require(arules)
require(arulesViz)


load("~/Projects/MBA/data/MASTER_MATRIX_PERMUTED.RData")


b.t.perm<-as(as.data.frame(m.final.perm),"transactions")
b.t<-b.t.perm
drug.idx<-grep("IC.50",labels(m.final.perm)[[2]])
drugs<-labels(m.final.perm)[[2]][drug.idx]
drugs.total<-c(paste(drugs,"Sensitive",sep="="),paste(drugs,"Resistant",sep="="))

rules<-apriori(b.t,parameter = list(support = (4/length(b.t)), confidence = 4/length(b.t), maxlen=2),appearance = list(rhs=drugs.total,default="lhs"))

temp.quality<-quality(rules)
temp.quality[,3]<-NA
temp.quality<-temp.quality[!duplicated(temp.quality), ]

for (i in 1:nrow(temp.quality)){
  print(i)
  rules.perm<-apriori(b.t.perm,parameter = list(support = temp.quality$support[i], confidence = temp.quality$confidence[i], maxlen=2))
  if (length(rules.perm)==0) next
  temp<-hist(log(quality(rules.perm)[,3]),breaks=80,plot=F)
  lift.estim.0.05<-exp(temp$breaks[which(cumsum(temp$counts)>sum(temp$counts)*0.95)[1]])
  temp.quality$lift[i]<-lift.estim.0.05
}
temp.quality$lift[temp.quality$lift<=1]<-NA
temp.quality$lift[is.na(temp.quality$lift)]<-min(temp.quality$lift,na.rm=T)


rm(b.t.perm)
rm(rules.perm)
gc()


status<-rep(NA,length(rules))
quality.rules<-quality(rules)
require(parallel)
require(arules)
cl<-makeCluster(4)
clusterExport(cl, "quality.rules")
clusterExport(cl, "temp.quality")
clusterExport(cl, "status")
for (i in 1:nrow(temp.quality)){
  print(i)
  clusterExport(cl, "i")
  temp.idx<-which(quality.rules$support==temp.quality$support[i] & quality.rules$confidence==temp.quality$confidence[i])
  clusterExport(cl, "temp.idx")
  status[temp.idx]<-parSapply (cl=cl, quality.rules$lift[temp.idx], function(lift) if(lift>temp.quality$lift[i]){temp.out<-"pass"}else{temp.out<-"fail"})
}
stopCluster(cl)
rm(cl)
gc()

rules.sig<-rules[which(status=="pass")]
fdr<-length(rules.sig)/length(rules)





#############################################
###############Train DataSet#################
#############################################

require(arules)
require(arulesViz)

load("~/Projects/MBA/data/Prediction/TRAIN.RData")
load("~/Projects/MBA/data/Prediction/TRAIN_PERMUTED.RData")


b.t<-as(as.data.frame(train),"transactions")
b.t.perm<-as(as.data.frame(train.perm),"transactions")
drug.idx<-grep("IC.50",labels(train)[[2]])
drugs<-labels(train)[[2]][drug.idx]
drugs.total<-c(paste(drugs,"Sensitive",sep="="),paste(drugs,"Resistant",sep="="))

rules<-apriori(b.t,parameter = list(support = (4/length(b.t)), confidence = 4/length(b.t), maxlen=2),appearance = list(rhs=drugs.total,default="lhs"))

temp.quality<-quality(rules)
temp.quality[,3]<-NA
temp.quality<-temp.quality[!duplicated(temp.quality), ]

for (i in 1:nrow(temp.quality)){
  print(i)
  rules.perm<-apriori(b.t.perm,parameter = list(support = temp.quality$support[i], confidence = temp.quality$confidence[i], maxlen=2))
  if (length(rules.perm)==0) next
  temp<-hist(log(quality(rules.perm)[,3]),breaks=80,plot=F)
  lift.estim.0.05<-exp(temp$breaks[which(cumsum(temp$counts)>sum(temp$counts)*0.95)[1]])
  temp.quality$lift[i]<-lift.estim.0.05
}
temp.quality$lift[temp.quality$lift<=1]<-NA
temp.quality$lift[is.na(temp.quality$lift)]<-min(temp.quality$lift,na.rm=T)


rm(b.t.perm)
rm(rules.perm)
gc()


status<-rep(NA,length(rules))
quality.rules<-quality(rules)
require(parallel)
require(arules)
cl<-makeCluster(4)
clusterExport(cl, "quality.rules")
clusterExport(cl, "temp.quality")
clusterExport(cl, "status")
for (i in 1:nrow(temp.quality)){
  print(i)
  clusterExport(cl, "i")
  temp.idx<-which(quality.rules$support==temp.quality$support[i] & quality.rules$confidence==temp.quality$confidence[i])
  clusterExport(cl, "temp.idx")
  status[temp.idx]<-parSapply (cl=cl, quality.rules$lift[temp.idx], function(lift) if(lift>temp.quality$lift[i]){temp.out<-"pass"}else{temp.out<-"fail"})
}
stopCluster(cl)
rm(cl)
gc()

rules.sig<-rules[which(status=="pass")]
save(rules.sig,file="~/Projects/MBA/data/Prediction/TRAIN_RULES_SIGNIFICANT_Sup4_Conf4in689_DYNAMIC_THRESHOLD_FDR0.05.RData")
write(rules.sig,file="~/Projects/MBA/data/Prediction/TRAIN_RULES_SIGNIFICANT_Sup4_Conf4in689_DYNAMIC_THRESHOLD_FDR0.05.tsv",sep='\t',row.names=F)
rm(rules)
gc()



#############################################
###############Drug2Drug#################
#############################################

require(arules)
require(arulesViz)

load("~/Projects/MBA/data/MASTER_MATRIX.RData")
load("~/Projects/MBA/data/MASTER_MATRIX_PERMUTED.RData")

b.t<-as(as.data.frame(m.final),"transactions")
b.t.perm<-as(as.data.frame(m.final.perm),"transactions")
drug.idx<-grep("IC.50",labels(m.final)[[2]])
drugs<-labels(m.final)[[2]][drug.idx]
drugs.total<-c(paste(drugs,"Sensitive",sep="="),paste(drugs,"Resistant",sep="="))

rules<-apriori(b.t,parameter = list(support = (4/length(b.t)), confidence = 4/length(b.t), maxlen=2))
rules.sub<- subset(rules, subset = (lhs %pin% "IC_50" & rhs %pin% "IC_50"))
rm(rules)
rules<-rules.sub
rm(rules.sub)

temp.quality<-quality(rules)
temp.quality[,3]<-NA
temp.quality<-temp.quality[!duplicated(temp.quality), ]
save(rules,file="~/Desktop/temp.RData")
rm(rules)
gc()

for (i in 1:nrow(temp.quality)){
  print(i)
  rules.perm<-apriori(b.t.perm,parameter = list(support = temp.quality$support[i], confidence = temp.quality$confidence[i], maxlen=2))
  if (length(rules.perm)==0) next
  temp<-hist(log(quality(rules.perm)[,3]),breaks=80,plot=F)
  lift.estim.0.05<-exp(temp$breaks[which(cumsum(temp$counts)>sum(temp$counts)*0.95)[1]])
  temp.quality$lift[i]<-lift.estim.0.05
}
temp.quality$lift[temp.quality$lift<=1]<-NA
temp.quality$lift[is.na(temp.quality$lift)]<-min(temp.quality$lift,na.rm=T)


rm(b.t.perm)
rm(rules.perm)
gc()

load("~/Desktop/temp.RData")
status<-rep(NA,length(rules))
quality.rules<-quality(rules)
require(arules)
for (i in 1:nrow(temp.quality)){
  print(i)
  temp.idx<-which(quality.rules$support==temp.quality$support[i] & quality.rules$confidence==temp.quality$confidence[i])
  status[temp.idx]<-sapply (quality.rules$lift[temp.idx], function(lift) if(lift>temp.quality$lift[i]){temp.out<-"pass"}else{temp.out<-"fail"})
}
gc()

rules.sig<-rules[which(status=="pass")]
save(rules.sig,file="~/Projects/MBA/data/RULES_SIGNIFICANT_DRUG2DRUG_Sup4_Conf4in689_DYNAMIC_THRESHOLD_FDR0.05.RData")
write(rules.sig,file="~/Projects/MBA/data/RULES_SIGNIFICANT_DRUG2DRUG_Sup4_Conf4in689_DYNAMIC_THRESHOLD_FDR0.05.tsv",sep='\t',row.names=F)
rm(rules)
gc()
#write.table(temp.quality,file="~/Projects/MBA/data/Thresholds_Sup_Conf_Lift_forFDR0.05_grid.tsv",sep='\t',row.names = F)
