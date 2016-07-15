########################################
###Create CCLP_Drug_Categorical.csv#####
########################################
cclp.drug<-read.csv("~/Projects/MBA/data/gdsc_manova_input_w5_proc.csv")
ic50.idx<-grep("_IC_50.1",names(cclp.drug))
b<-cclp.drug[,ic50.idx]
b2<-b
b2<-apply(b2,2,scale,scale=T,center=T)
b3<-b2
b3[b2>=1]<-"Resistant"
b3[b2<=(-1)]<-"Sensitive"
b3[b2<1 & b2>-1]<-NA
b3<-cbind(Cell.Line=as.character(cclp.drug[,1]),b3)
write.table(b3,file="~/Projects/MBA/data/CCLP_Drug_Categorical.tsv",sep='\t',row.names = F)



########################################
###Create MASTER MATRIX#####
########################################

##CCLP_Mutation
a<-read.csv("~/Projects/MBA/data/CosmicCLP_MutantExport.tsv",sep='\t')
a$Sample.name<-gsub(" ","",a$Sample.name)
a$Sample.name<-gsub("-","",a$Sample.name)
a$Sample.name<-gsub("\\.","",a$Sample.name)
a$Sample.name<-gsub("/","",a$Sample.name)
a$Sample.name<-toupper(a$Sample.name)
a$Gene.name<-as.character(gsub("_.*","",a$Gene.name))
genes.mut<-unique(a$Gene.name)

##CCLP_Drug
cclp.drug<-read.csv("~/Projects/MBA/data/CCLP_Drug_Categorical.tsv",sep='\t')
cclp.drug[,1]<-gsub(" ","",cclp.drug[,1])
cclp.drug[,1]<-gsub("-","",cclp.drug[,1])
cclp.drug[,1]<-gsub("\\.","",cclp.drug[,1])
cclp.drug[,1]<-gsub("/","",cclp.drug[,1])
cclp.drug[,1]<-toupper(cclp.drug[,1])
cell.lines<-unique(cclp.drug[,1])
compounds<-names(cclp.drug)[-1]

##CCLP_GE
cclp.ge<-read.csv("~/Projects/MBA/data/CCLP_CompleteGeneExpression.tsv",sep='\t')
cclp.ge$SAMPLE_NAME<-gsub(" ","",cclp.ge$SAMPLE_NAME)
cclp.ge$SAMPLE_NAME<-gsub("-","",cclp.ge$SAMPLE_NAME)
cclp.ge$SAMPLE_NAME<-gsub("\\.","",cclp.ge$SAMPLE_NAME)
cclp.ge$SAMPLE_NAME<-gsub("/","",cclp.ge$SAMPLE_NAME)
cclp.ge$SAMPLE_NAME<-toupper(cclp.ge$SAMPLE_NAME)
genes.ge<-unique(as.character(cclp.ge$GENE_NAME))

#CCLP_CNV
cclp.cnv<-read.csv("~/Projects/MBA/data/CCLP_CNV.tsv",sep='\t')
cclp.cnv$SAMPLE_NAME<-gsub(" ","",cclp.cnv$SAMPLE_NAME)
cclp.cnv$SAMPLE_NAME<-gsub("-","",cclp.cnv$SAMPLE_NAME)
cclp.cnv$SAMPLE_NAME<-gsub("\\.","",cclp.cnv$SAMPLE_NAME)
cclp.cnv$SAMPLE_NAME<-gsub("/","",cclp.cnv$SAMPLE_NAME)
cclp.cnv$SAMPLE_NAME<-toupper(cclp.cnv$SAMPLE_NAME)
genes.cnv<-unique(as.character(cclp.cnv$GENE_NAME))

m.mut<-matrix(NA,nrow=length(cell.lines),ncol=length(genes.mut),dimnames=list(cell.lines,genes.mut))
m.drug<-matrix(NA,nrow=length(cell.lines),ncol=length(compounds),dimnames=list(cell.lines,compounds))
m.ge<-matrix(NA,nrow=length(cell.lines),ncol=length(genes.ge),dimnames=list(cell.lines,genes.ge))
m.cnv<-matrix(NA,nrow=length(cell.lines),ncol=length(genes.cnv),dimnames=list(cell.lines,genes.cnv))

for (i in 1:nrow(a)){
  if (length(which(cell.lines==a$Sample.name[i]))>0){
    if (a$Mutation.Description[i]!="Substitution - coding silent" & is.na(m.mut[as.character(a$Sample.name[i]),as.character(a$Gene.name[i])])){
      m.mut[as.character(a$Sample.name[i]),as.character(a$Gene.name[i])]<-"Mut"  
    } else {next}
    
    }
}

m.drug<-as.matrix(cclp.drug[,-1])

for (i in 1:nrow(cclp.ge)){
  if (length(which(cell.lines==cclp.ge$SAMPLE_NAME[i]))>0) m.ge[as.character(cclp.ge$SAMPLE_NAME[i]),as.character(cclp.ge$GENE_NAME[i])]<-as.character(cclp.ge$REGULATION[i])
}
m.ge[m.ge=='normal']<-NA


for (i in 1:nrow(cclp.cnv)){
  if (length(which(cell.lines==cclp.cnv$SAMPLE_NAME[i]))>0) {
    if (is.na(m.cnv[as.character(cclp.cnv$SAMPLE_NAME[i]),as.character(cclp.cnv$GENE_NAME[i])])){
      m.cnv[as.character(cclp.cnv$SAMPLE_NAME[i]),as.character(cclp.cnv$GENE_NAME[i])]<-as.character(cclp.cnv$MUT_TYPE[i])  
    } else if (m.cnv[as.character(cclp.cnv$SAMPLE_NAME[i]),as.character(cclp.cnv$GENE_NAME[i])]=="GAIN" & as.character(cclp.cnv$MUT_TYPE[i])=="GAIN") {
    next  
    } else if (m.cnv[as.character(cclp.cnv$SAMPLE_NAME[i]),as.character(cclp.cnv$GENE_NAME[i])]=="LOSS" & as.character(cclp.cnv$MUT_TYPE[i])=="LOSS") {
    next  
    } else if ((m.cnv[as.character(cclp.cnv$SAMPLE_NAME[i]),as.character(cclp.cnv$GENE_NAME[i])]=="GAIN" & as.character(cclp.cnv$MUT_TYPE[i])=="LOSS")){
      m.cnv[as.character(cclp.cnv$SAMPLE_NAME[i]),as.character(cclp.cnv$GENE_NAME[i])]<-"GAIN,LOSS"  
    } else if ((m.cnv[as.character(cclp.cnv$SAMPLE_NAME[i]),as.character(cclp.cnv$GENE_NAME[i])]=="LOSS" & as.character(cclp.cnv$MUT_TYPE[i])=="GAIN")){
      m.cnv[as.character(cclp.cnv$SAMPLE_NAME[i]),as.character(cclp.cnv$GENE_NAME[i])]<-"GAIN,LOSS"
    }
    
   }
}



m.total<-cbind(m.mut,m.cnv,m.ge,m.drug)
dimnames(m.total)<-list(cell.lines,c(genes.mut,paste(genes.cnv,"CNV",sep="_"),paste(genes.ge,"GE",sep="_"),compounds))


cell.lines.to.remove<-c("DMS153","CALU1","BOKU","DV90","NCIH1618","SCLC21H","ECC4","K052","JRT3T35","LC1F","MDAMB134VI","NTERASCLD1","NCIH1522","TI73","EFE184","NCIH719","NCIH2107","PLCPRF5")
rm.idx<-match(cell.lines.to.remove,labels(m.total)[[1]])
m.final<-m.total[-rm.idx,]
save(m.final,file="~/Projects/MBA/data/MASTER_MATRIX.RData")
load("~/Projects/MBA/data/MASTER_MATRIX.RData")

##Add Tissue info to m.final from cclp.drug
Tissue<-NA
cclp.drug$Cell.Line<-gsub(" ","",cclp.drug$Cell.Line)
cclp.drug$Cell.Line<-gsub("-","",cclp.drug$Cell.Line)
cclp.drug$Cell.Line<-gsub("\\.","",cclp.drug$Cell.Line)
cclp.drug$Cell.Line<-gsub("/","",cclp.drug$Cell.Line)
cclp.drug$Cell.Line<-toupper(cclp.drug$Cell.Line)
for (i in 1:nrow(m.final)){
  idx<-which(cclp.drug$Cell.Line==labels(m.final)[[1]][i])
  Tissue[i]<-as.character(cclp.drug$Tissue[idx])
}
m.final<-cbind(Tissue,m.final)
save(m.final,file="~/Projects/MBA/data/MASTER_MATRIX.RData")



###Do permutations (Parallel Computing)
require(parallel)
cl<-makeCluster(7)
clusterExport(cl, "m.final")
m.final.perm<-parSapply (cl=cl, 1:ncol(m.final), function (col) m.final[,col]<-sample(m.final[,col]))
dimnames(m.final.perm)<-list(labels(m.final)[[1]],labels(m.final)[[2]])
stopCluster(cl)
save(m.final.perm,file="~/Projects/MBA/data/MASTER_MATRIX_PERMUTED.RData")







###Make Train & Test, normal & permuted ###
load("~/Projects/MBA/data/MASTER_MATRIX.RData")
load("~/Projects/MBA/data/MASTER_MATRIX_PERMUTED.RData")

tissue.table<-sort(table(m.final[,1]),decreasing=T)
m.final.row.names<-labels(m.final)[[1]]
m.final.perm.row.names<-labels(m.final.perm)[[1]]
for (i in 1:length(tissue.table)){
  temp.idx<-which(m.final[,1]==names(tissue.table[i])) 
  temp.idx.perm<-which(m.final.perm[,1]==names(tissue.table[i])) 
  temp.m.final<-m.final[temp.idx,]
  temp.m.final.perm<-m.final.perm[temp.idx.perm,]
  temp.m.final.row.names<-m.final.row.names[which(m.final[,1]==names(tissue.table[i])) ]
  temp.m.final.perm.row.names<-m.final.perm.row.names[which(m.final.perm[,1]==names(tissue.table[i])) ]
  train.no<-round((tissue.table[i]/3)*2,digits = 0)
  train.random.idx<-sample(c(1:tissue.table[i]),train.no)
  print(paste(tissue.table[i],train.no,paste(train.random.idx,collapse="-"),sep=","))
  if (i==1){
    train<-temp.m.final[train.random.idx,]
    test<-temp.m.final[-train.random.idx,]
    train.perm<-temp.m.final.perm[train.random.idx,]
    test.perm<-temp.m.final.perm[-train.random.idx,]
    train.row.names<-temp.m.final.row.names[train.random.idx]
    test.row.names<-temp.m.final.row.names[-train.random.idx]
    train.perm.row.names<-temp.m.final.perm.row.names[train.random.idx]
    test.perm.row.names<-temp.m.final.perm.row.names[-train.random.idx]
    next
  }
  if (tissue.table[i]==1){
    train<-rbind(train,temp.m.final)
    train.perm<-rbind(train.perm,temp.m.final.perm)
    train.row.names<-c(train.row.names,temp.m.final.row.names)
    train.perm.row.names<-c(train.perm.row.names,temp.m.final.perm.row.names)
    next
  } else {
    train.no<-round((tissue.table[i]/3)*2,digits = 0)
    train.random.idx<-sample(c(1:tissue.table[i]),train.no)
    train<-rbind(train,temp.m.final[train.random.idx,])
    test<-rbind(test,temp.m.final[-train.random.idx,])
    train.perm<-rbind(train.perm,temp.m.final.perm[train.random.idx,])
    test.perm<-rbind(test.perm,temp.m.final.perm[-train.random.idx,])
    train.row.names<-c(train.row.names,temp.m.final.row.names[train.random.idx])
    test.row.names<-c(test.row.names,temp.m.final.row.names[-train.random.idx])
    train.perm.row.names<-c(train.perm.row.names,temp.m.final.perm.row.names[train.random.idx])
    test.perm.row.names<-c(test.perm.row.names,temp.m.final.perm.row.names[-train.random.idx])
  }
}

dimnames(train)<-list(train.row.names,labels(m.final)[[2]])
dimnames(test)<-list(test.row.names,labels(m.final)[[2]])
dimnames(train.perm)<-list(train.perm.row.names,labels(m.final.perm)[[2]])
dimnames(test.perm)<-list(test.perm.row.names,labels(m.final.perm)[[2]])

save(train,file="~/Projects/MBA/data/Prediction/TRAIN.RData")
save(test,file="~/Projects/MBA/data/Prediction/TEST.RData")
save(train.perm,file="~/Projects/MBA/data/Prediction/TRAIN_PERMUTED.RData")
save(test.perm,file="~/Projects/MBA/data/Prediction/TEST_PERMUTED.RData")


train.table<-sort(table(train[,1]),decreasing=T)
test.table<-sort(table(test[,1]),decreasing=T)

train.perm.table<-sort(table(train.perm[,1]),decreasing=T)
test.perm.table<-sort(table(test.perm[,1]),decreasing=T)

#######################################################################################
####Produce Train & Test matrices were Gene-expression factor levels have been#########
####replaced by their original z-transformed gene expression values
#######################################################################################
cclp.ge<-read.csv("~/Projects/MBA/data/CCLP_CompleteGeneExpression.tsv",sep='\t')
cclp.ge$SAMPLE_NAME<-gsub(" ","",cclp.ge$SAMPLE_NAME)
cclp.ge$SAMPLE_NAME<-gsub("-","",cclp.ge$SAMPLE_NAME)
cclp.ge$SAMPLE_NAME<-gsub("\\.","",cclp.ge$SAMPLE_NAME)
cclp.ge$SAMPLE_NAME<-gsub("/","",cclp.ge$SAMPLE_NAME)
cclp.ge$SAMPLE_NAME<-toupper(cclp.ge$SAMPLE_NAME)
for (i in 1:nrow(cclp.ge)){
  print(i)
  try(train[cclp.ge$SAMPLE_NAME[i],paste(cclp.ge$GENE_NAME[i],"GE",sep="_")]<-cclp.ge$Z_SCORE[i],silent=T)
  try(test[cclp.ge$SAMPLE_NAME[i],paste(cclp.ge$GENE_NAME[i],"GE",sep="_")]<-cclp.ge$Z_SCORE[i],silent=T)
}
save(train,file="~/Projects/MBA/data/Prediction/TRAIN_GE_NUM.RData")
save(test,file="~/Projects/MBA/data/Prediction/TEST_GE_NUM.RData")



#######################################################################################
####Produce Train & Test matrices were Gene-expression factor levels have been#########
####replaced by their original z-transformed gene expression values and IC50 values
#######################################################################################
load("~/Projects/MBA/data/Prediction/TRAIN_GE_NUM.RData")
load("~/Projects/MBA/data/Prediction/TEST_GE_NUM.RData")
cclp.drug<-read.csv("~/Projects/MBA/data/gdsc_manova_input_w5_proc.csv")
ic50.idx<-grep("_IC_50.1",names(cclp.drug))
b<-cclp.drug[,ic50.idx]
b2<-b
b2<-apply(b2,2,scale,scale=T,center=T)
dimnames(b2)<-list(toupper(gsub("-","",cclp.drug[,1])),dimnames(b2)[[2]])
for (col in labels(b2)[[2]]){
train[,col]<-as.numeric(train[,col])
for (row in labels(train)[[1]]){
train[row,col]<-b2[row,col]  
}  
}
for (col in labels(b2)[[2]]){
  test[,col]<-as.numeric(test[,col])
  for (row in labels(test)[[1]]){
    test[row,col]<-b2[row,col]  
  }  
}

######################################################
##Make drug matirx with actual z-transformed IC50 values
#########################################################
cclp.drug<-read.csv("~/Projects/MBA/data/gdsc_manova_input_w5_proc.csv")
ic50.idx<-grep("_IC_50.1",names(cclp.drug))
b<-cclp.drug[,ic50.idx]
b2<-b
b3<-apply(b2,2,scale,scale=T,center=T)
#b3<-b2
Cell.Line=as.character(cclp.drug[,1])
Cell.Line<-gsub(" ","",Cell.Line)
Cell.Line<-gsub("-","",Cell.Line)
Cell.Line<-gsub("\\.","",Cell.Line)
Cell.Line<-gsub("/","",Cell.Line)
Cell.Line<-toupper(Cell.Line)
dimnames(b3)<-list(Cell.Line,labels(b3)[[2]])


load("~/Projects/MBA/data/MASTER_MATRIX.RData")
m.idx<-match(labels(b3)[[1]],labels(m.final.drugs)[[1]])
m.final.drugs<-b3[-which(is.na(m.idx)),]
save(m.final.drugs,file="~/Projects/MBA/data/Matrix_Drugs_IC50.RData")





###Produce processed version of rules###
########################################




##One-way##
rules4<-read.csv("~/Projects/MBA/data/RULES_SIGNIFICANT_Sup4_Conf4in689_DYNAMIC_THRESHOLD_FDR0.05.tsv",sep='\t')
rules4[,1]<-as.character(rules4[,1])
out4<-data.frame(lhs=rep(NA,nrow(rules4)),
                 rhs=rep(NA,nrow(rules4)),
                 support=rep(NA,nrow(rules4)),
                 confidence=rep(NA,nrow(rules4)),
                 lift=rep(NA,nrow(rules4))
)
rules4[,1]<-gsub("\\{|\\}","",rules4[,1])
tmp.split<-strsplit(rules4[,1]," => ")
for (i in 1:nrow(out4)){
  out4[i,1:2]<-unlist(tmp.split[i])
  out4[i,3:5]<-rules4[i,2:4]
}
save(out4,file="~/Projects/MBA/data/rules4_proc.RData")



##Two way
rules6.2way<-read.csv("~/Projects/MBA/data/TRAIN_RULES_TWO_WAY_SIGNIFICANT_Sup6_Conf6in689_DYNAMIC_THRESHOLD_FDR0.05.tsv",sep='\t')
rules6.2way[,1]<-as.character(rules6.2way[,1])
out6.2way<-data.frame(lhs=rep(NA,nrow(rules6.2way)),
                      rhs=rep(NA,nrow(rules6.2way)),
                      support=rep(NA,nrow(rules6.2way)),
                      confidence=rep(NA,nrow(rules6.2way)),
                      lift=rep(NA,nrow(rules6.2way))
)
rules6.2way[,1]<-gsub("\\{|\\}","",rules6.2way[,1])
tmp.split<-strsplit(rules6.2way[,1]," => ")
out6.2way[,1:2]<-t(sapply(tmp.split,unlist))
out6.2way[,3:5]<-rules6.2way[,2:4]
save(out6.2way,file="~/Projects/MBA/data/train_rules6_2way_proc.RData")



###############
##########Train
###############
rules4<-read.csv("~/Projects/MBA/data/Prediction/TRAIN_RULES_SIGNIFICANT_Sup4_Conf4in689_DYNAMIC_THRESHOLD_FDR0.05.tsv",sep='\t')
rules4[,1]<-as.character(rules4[,1])
out4.train<-data.frame(lhs=rep(NA,nrow(rules4)),
                       rhs=rep(NA,nrow(rules4)),
                       support=rep(NA,nrow(rules4)),
                       confidence=rep(NA,nrow(rules4)),
                       lift=rep(NA,nrow(rules4))
)
rules4[,1]<-gsub("\\{|\\}","",rules4[,1])
tmp.split<-strsplit(rules4[,1]," => ")
for (i in 1:nrow(out4.train)){
  out4.train[i,1:2]<-unlist(tmp.split[i])
  out4.train[i,3:5]<-rules4[i,2:4]
}
save(out4.train,file="~/Projects/MBA/data/Prediction/rules4_train_proc.RData")
