par(family="Times") ##SOOOOOS TIMES NEW ROMAN
require(arules)
require(arulesViz)

##CCLP_Drug
cclp.drug<-read.csv("~/Projects/MBA/data/CCLP_Drug_Categorical.tsv",sep='\t')
cclp.drug[,1]<-gsub(" ","",cclp.drug[,1])
cclp.drug[,1]<-gsub("-","",cclp.drug[,1])
cclp.drug[,1]<-gsub("\\.","",cclp.drug[,1])
cclp.drug[,1]<-gsub("/","",cclp.drug[,1])
cclp.drug[,1]<-toupper(cclp.drug[,1])


load("~/Projects/MBA/data/RULES_SIGNIFICANT_Sup4_Conf4in689_DYNAMIC_THRESHOLD_FDR0.05.RData")
#load("/home/kvougas/Projects/MBA/data/RULES_2WAY_SIGNIFICANT_SUP7_CONF01_DYNAMIC_THRESHOLD_FDR0.05.RData")
load("~/Projects/MBA/data/MASTER_MATRIX.RData")






##All
rules.sig.res<- subset(rules.sig, subset = (rhs %pin% "Resistant"))
rules.sig.sens<- subset(rules.sig, subset = (rhs %pin% "Sensi"))
rules.sig.sens<-sort(rules.sig.sens,by="support")
#rules.sig.sens<-sort(rules.sig.sens,by="lift")
#rules.sig.sens<-sort(rules.sig.sens,by="confidence")
plot(rules.sig.sens[1:2000], method="grouped", interactive = F,control=list(k=200))

rules.sig.res<-sort(rules.sig.res,by="support")
#rules.sig.res<-sort(rules.sig.res,by="lift")
#rules.sig.res<-sort(rules.sig.res,by="confidence")
plot(rules.sig.res[1:2000], method="grouped", interactive = F,control=list(k=200))
##


##Tissue
rules.sig.res<- sort(subset(rules.sig, subset = (lhs %pin% "Tissue=" & rhs %pin% "Resistant" & confidence>0.5)),by="confidence")
rules.sig.sens<- sort(subset(rules.sig, subset = (lhs %pin% "Tissue=" & rhs %pin% "Sensi"& confidence>0.5)),by="confidence")
rules.sig.sens<-sort(rules.sig.sens,by="support")
#rules.sig.sens<-sort(rules.sig.sens,by="lift")
#rules.sig.sens<-sort(rules.sig.sens,by="confidence")
plot(rules.sig.sens, method="grouped", interactive = F,control=list(k=50))

rules.sig.res<-sort(rules.sig.res,by="support")
#rules.sig.res<-sort(rules.sig.res,by="lift")
#rules.sig.res<-sort(rules.sig.res,by="confidence")
plot(rules.sig.res, method="grouped", interactive = F,control=list(k=50))
##




#BRAF-MEK-PI3K
idx<-scan("~/Projects/MBA/data/BRAF_MEK_CORRECT_SENS_RULES4_IDX.txt")
#idx<-scan("~/Projects/MBA/data/PI3K_SENS_RULES4_IDX.txt")
rules.sub<-rules.sig[idx] #Select either of the two from lines 56,57
rules.sub<-sort(rules.sub,by="support")
#rules.sub<-sort(rules.sub,by="lift")
#rules.sub<-sort(rules.sub,by="confidence")
temp<-plot(rules.sub, method="grouped", interactive = F,control=list(k=50))



#IC50-Heatmaps
require(arules)
require(gplots)
load("~/Projects/MBA/data/RULES_SIGNIFICANT_Sup4_Conf4in689_DYNAMIC_THRESHOLD_FDR0.05.RData")
load("~/Projects/MBA/data/MASTER_MATRIX.RData")
load("~/Projects/MBA/data/Matrix_Drugs_IC50.RData")



#BRAF-MEK
idx<-scan("~/Projects/MBA/data/BRAF_MEK_SENS_RULES4_IDX.txt")
rules.sub<-rules.sig[idx]
rules.sub<- subset(rules.sub, subset = (lhs %pin% "BRAF=Mut" & rhs %pin% "Sensitive"))
temp.drugs<-gsub("\\{|\\}|=Sensitive","",as.character(unlist(inspect(rhs(rules.sub)))))
idx<-which(m.final[,'BRAF']=="Mut")
temp.m.final.drugs<-m.final.drugs[idx,]
m.idx<-match(temp.drugs,labels(temp.m.final.drugs)[[2]])
temp.m.final.drugs<-temp.m.final.drugs[,m.idx]
temp.m.final.drugs[is.na(temp.m.final.drugs)]<-0
dimnames(temp.m.final.drugs)<-list(labels(temp.m.final.drugs)[[1]],gsub("_IC_50\\.1","",labels(temp.m.final.drugs)[[2]]))
jpeg("~/Projects/MBA/data/Viz/BRAF_Mut_Heatmap.jpg",width=2400,height=1820,pointsize = 12)
heatmap.2(temp.m.final.drugs,col=bluered(100),symm=F,symkey=F,symbreaks=T, scale="none",cexCol=3,cexRow=1,trace="none",key=F,srtCol = 10,srtRow = 33,dendrogram="none")
dev.off()

#PI3K
idx<-scan("~/Projects/MBA/data/PI3K_SENS_RULES4_IDX.txt")
rules.sub<-rules.sig[idx]
rules.sub<- subset(rules.sub, subset = (lhs %pin% "PIK3CA=Mut" & rhs %pin% "Sensitive"))
temp.drugs<-gsub("\\{|\\}|=Sensitive","",as.character(unlist(inspect(rhs(rules.sub)))))
idx<-which(m.final[,'PIK3CA']=="Mut")
temp.m.final.drugs<-m.final.drugs[idx,]
m.idx<-match(temp.drugs,labels(temp.m.final.drugs)[[2]])
temp.m.final.drugs<-temp.m.final.drugs[,m.idx]
temp.m.final.drugs[is.na(temp.m.final.drugs)]<-0
dimnames(temp.m.final.drugs)<-list(labels(temp.m.final.drugs)[[1]],gsub("_IC_50\\.1","",labels(temp.m.final.drugs)[[2]]))
jpeg("~/Projects/MBA/data/Viz/PIK3CA_Mut_Heatmap.jpg",width=2400,height=1820,pointsize = 12)
heatmap.2(temp.m.final.drugs,col=bluered(100),symm=F,symkey=F,symbreaks=T, scale="none",cexCol=3,cexRow=1,trace="none",key=F,srtCol = 10,srtRow = 33,dendrogram="none")
dev.off()

#Tissue=lung_small_cell_carcinoma
load("~/Projects/MBA/data/RULES_SIGNIFICANT_Sup4_Conf4in689_DYNAMIC_THRESHOLD_FDR0.05.RData")
rules.sub<- subset(rules.sig, subset = (lhs %pin% "Tissue=lung_small_cell_carcinoma" & rhs %pin% "Resistant"))
#write(rules.sub,file="/home/kvougas/Projects/MBA/data/FMO5_Res_1way.tsv",sep='\t')
temp.drugs<-gsub("\\{|\\}|=Resistant","",as.character(unlist(inspect(rhs(rules.sub)))))
idx<-which(m.final[,'Tissue']=="lung_small_cell_carcinoma")
temp.m.final.drugs<-m.final.drugs[idx,]
m.idx<-match(temp.drugs,labels(temp.m.final.drugs)[[2]])
temp.m.final.drugs<-temp.m.final.drugs[,m.idx]
temp.m.final.drugs[is.na(temp.m.final.drugs)]<-0
dimnames(temp.m.final.drugs)<-list(labels(temp.m.final.drugs)[[1]],gsub("_IC_50\\.1","",labels(temp.m.final.drugs)[[2]]))
jpeg("~/Projects/MBA/data/Viz/Tissue_Small_Cell_Lung_Carcinoma.jpg",width=2400,height=1820,pointsize = 15)
heatmap.2(temp.m.final.drugs,col=bluered(100),symm=F,symkey=F,symbreaks=T, scale="none",cexCol=0.9,cexRow=0.9,trace="none",key=F,srtCol = 20,srtRow = 33,dendrogram="none")
#heatmap.2(temp.m.final.drugs,col=bluered(100),symm=F,symkey=F,symbreaks=T, scale="none",cexCol=3,cexRow=1,trace="none",key=F,srtCol = 10,srtRow = 33,dendrogram="none")
dev.off()
#




##Boxplots for prediction performance
require(ggplot2)
dl.rf.comb<-read.csv("~/Projects/MBA/data/DL_RF_Comb_for_ggplot2.txt",sep='\t')
#Sens
p <- ggplot(dl.rf.comb, aes(x=Classifier, y=Sensitivity)) 
 p + 
  geom_boxplot(fill = c("#0000ff", "#8b0000"),notch=F)+ 
  geom_jitter(shape=16, position=position_jitter(0.2))+ 
  theme(panel.background = element_rect(fill='white'))

#ACC
 p <- ggplot(dl.rf.comb, aes(x=Classifier, y=ACC)) 
 p + 
   geom_boxplot(fill = c("#0000ff", "#8b0000"),notch=F)+ 
   geom_jitter(shape=16, position=position_jitter(0.2))+ 
   theme(panel.background = element_rect(fill='white'))
 
 #NPV
 p <- ggplot(dl.rf.comb, aes(x=Classifier, y=NPV)) 
 p + 
   geom_boxplot(fill = c("#0000ff", "#8b0000"),notch=F)+ 
   geom_jitter(shape=16, position=position_jitter(0.2))+ 
   theme(panel.background = element_rect(fill='white'))
 
 #Specificity
 p <- ggplot(dl.rf.comb, aes(x=Classifier, y=Specificity)) 
 p + 
   geom_boxplot(fill = c("#0000ff", "#8b0000"),notch=F)+ 
   geom_jitter(shape=16, position=position_jitter(0.2))+ 
   theme(panel.background = element_rect(fill='white'))
 
 #PPV
 p <- ggplot(dl.rf.comb, aes(x=Classifier, y=PPV)) 
 p + 
   geom_boxplot(fill = c("#0000ff", "#8b0000"),notch=F)+ 
   geom_jitter(shape=16, position=position_jitter(0.2))+ 
   theme(panel.background = element_rect(fill='white'))
 
 #FPR
 p <- ggplot(dl.rf.comb, aes(x=Classifier, y=FPR)) 
 p + 
   geom_boxplot(fill = c("#0000ff", "#8b0000"),notch=F)+ 
   geom_jitter(shape=16, position=position_jitter(0.2))+ 
   theme(panel.background = element_rect(fill='white'))
 
