###############################################################
################START FROM HERE
###############################################################

par(family="Times") ##SOOOOOS TIMES NEW ROMAN
load("~/Projects/MBA/data/Prediction/rules4_train_proc.RData")
load("~/Projects/MBA/data/Prediction/TRAIN_GE_NUM.RData")
load("~/Projects/MBA/data/Prediction/TEST_GE_NUM.RData")
train.ge.idx<-grep("_GE",labels(train)[[2]])
test.ge.idx<-grep("_GE",labels(test)[[2]])
temp.train<-train[,train.ge.idx]
temp.test<-test[,test.ge.idx]
temp.train[is.na(temp.train)]<-"0"
temp.test[is.na(temp.test)]<-"0"
train[,train.ge.idx]<-temp.train
test[,test.ge.idx]<-temp.test
temp.test<-test[,2:39615] 
temp.train<-train[,2:39615] 
temp.test[is.na(temp.test)]<-"Normal"
temp.train[is.na(temp.train)]<-"Normal"
test[,2:39615]<-temp.test 
train[,2:39615]<-temp.train 
train<-as.data.frame(train)
test<-as.data.frame(test)
train.ge.idx<-grep("_GE",names(train))
test.ge.idx<-grep("_GE",names(test))
for (i in train.ge.idx){
  train[,i]<-as.numeric(levels(train[,i]))[train[,i]]
}
for (i in test.ge.idx){
  test[,i]<-as.numeric(levels(test[,i]))[test[,i]]
}
#Remove rows which do not contain gene expression info
train.rm.idx<-which(train[,28000]==0)
test.rm.idx<-which(test[,28000]==0)
train<-train[-train.rm.idx,]
test<-test[-test.rm.idx,]

##################
##################




##Prediction and evalutation
require(ROCR)

require(h2o)
localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, nthreads=5, max_mem_size="9G")


test.genes.only<-test[,2:39615]
drugs<-unique(out4.train$rhs)

result<-data.frame(Drug=rep(NA,length(drugs)),
                   Drug_Status=rep(NA,length(drugs)),
                   Predictor_Genes=rep(NA,length(drugs)),
                   Predictor_Genes_No=rep(NA,length(drugs))
)

counter=1
for (drug in drugs){
  
  temp.rules<-out4.train[out4.train$rhs==drug,]
  
  temp.rules2<-temp.rules[temp.rules$support>=quantile(temp.rules$support)[2],]#Remove rows which do not contain gene expression info
  
  temp.rules2$lhs<-gsub("=under|=over|=GAIN|=LOSS|=Mut","",temp.rules2$lhs)
  predictor.genes.0<-unique(unlist(strsplit(temp.rules2$lhs,",")))
  if (length(grep("=",predictor.genes.0))>0) predictor.genes.0<-predictor.genes.0[-grep("=",predictor.genes.0)]
  if (length(predictor.genes.0)==0) next
  no.of.rounds<-round(length(predictor.genes.0)/200,digits=0)*3
  if (no.of.rounds<1) no.of.rounds<-3  
  print(paste(drug,no.of.rounds,length(predictor.genes.0),sep=","))
  temp.auc<-NA
  temp.auc.counter<-1
  for (bagging in 1:no.of.rounds){
    if (length(predictor.genes.0)>=200) {
      predictor.genes<-sample(predictor.genes.0,200)
    } else {
      predictor.genes<-sample(predictor.genes.0,round(length(predictor.genes.0)*0.75,digits=0))
    }
    drug2<-unlist(strsplit(drug,"="))[1]
    drug.status<-unlist(strsplit(drug,"="))[2]
    temp.train<-train[,c(drug2,predictor.genes,"Tissue")]
    temp.test<-test[,c(drug2,predictor.genes,"Tissue")]
    
    if(drug.status=="Sensitive"){
      temp.train[,1]<-as.character(temp.train[,1])
      temp.test[,1]<-as.character(temp.test[,1])
      temp.train[is.na(temp.train[,1]),1]<-"Non-Sensitive"
      temp.train[temp.train[,1]=="Resistant",1]<-"Non-Sensitive"
      temp.test[is.na(temp.test[,1]),1]<-"Non-Sensitive"
      temp.test[temp.test[,1]=="Resistant",1]<-"Non-Sensitive"
      temp.train[,1]<-as.factor(temp.train[,1])
      temp.test[,1]<-as.factor(temp.test[,1])
    } else {
      temp.train[,1]<-as.character(temp.train[,1])
      temp.test[,1]<-as.character(temp.test[,1])
      temp.train[is.na(temp.train[,1]),1]<-"Non-Resistant"
      temp.train[temp.train[,1]=="Sensitive",1]<-"Non-Resistant"
      temp.test[is.na(temp.test[,1]),1]<-"Non-Resistant"
      temp.test[temp.test[,1]=="Sensitive",1]<-"Non-Resistant"
      temp.train[,1]<-as.factor(temp.train[,1])
      temp.test[,1]<-as.factor(temp.test[,1])
    }
    
    ge.idx<-grep("_GE",labels(temp.train)[[2]])
    right.columns<-c(1:ncol(temp.train))
    right.columns<-right.columns[-ge.idx]
    right.columns<-right.columns[-1]
    right.columns<-right.columns[-length(right.columns)]
    
    for(i in right.columns){
      temp.train[,i]<-as.character(temp.train[,i])
      temp.test[,i]<-as.character(temp.test[,i])
      train.levels<-sort(unique(temp.train[,i]))
      test.levels<-sort(unique(temp.test[,i]))
      total.levels<-sort(unique(c(train.levels,test.levels)))
      m.idx.train<-match(total.levels,train.levels)
      missing.train<-total.levels[is.na(m.idx.train)]
      if (length(missing.train)>0){
        for (impute.value in missing.train){
          normal.idx<-which(temp.train[,i]=="Normal")
          temp.train[sample(normal.idx,1),i]<-impute.value
          #temp.train[sample(c(1:nrow(temp.train)),1),i]<-impute.value
        }
      }
      m.idx.test<-match(total.levels,test.levels)
      missing.test<-total.levels[is.na(m.idx.test)]
      if (length(missing.test)>0){
        for (impute.value in missing.test){
          normal.idx<-which(temp.test[,i]=="Normal")
          temp.test[sample(normal.idx,1),i]<-impute.value
          #temp.test[sample(c(1:nrow(temp.test)),1),i]<-impute.value
        }
      }
      temp.train[,i]<-as.factor(temp.train[,i])
      temp.test[,i]<-as.factor(temp.test[,i])  
    }
    
    
    temp.train.h2o<-as.h2o(temp.train,destination_frame = "temp.train.h2o")
    temp.test.h2o<-as.h2o(temp.test,destination_frame = "temp.test.h2o")
    
    
    model1 <- 
      h2o.randomForest(x = 2:ncol(temp.train),  
                       y = 1,   
                       training_frame = temp.train.h2o,
                       ntrees=round(length(predictor.genes)/2,digits=0),
                       balance_classes=T,
                       nfolds = 5 
                      ) 
    temp.auc[temp.auc.counter]<-model1@model$cross_validation_metrics@metrics$AUC
    temp.auc.counter=temp.auc.counter+1
    
    p.h2o.1 <- h2o.predict(model1, temp.test.h2o)
    p1<-as.data.frame(p.h2o.1)
    p1<-p1[,-1]
    p<-p1
    if (bagging==1) {p.comb<-p} else {p.comb<-cbind(p.comb,p)}
  }
  p.final<-p.comb[,seq(2,ncol(p.comb),2)]
  p.final<-apply(p.final,1,weighted.mean,temp.auc)
  
  
  if (sum(p.final,na.rm=T)==0) next
  
  pred <- prediction(p.final, temp.test[,1])
  perf <- performance(pred,"tpr","fpr")
  
  print("CHECK")
  filename<-paste("~/Projects/MBA/data/Prediction_h2o_rf_quantile_support_filter_GE_NUM/plots/",drug,".jpg",sep="")
  jpeg(filename,width=640,height=640,pointsize=12)
  plot(perf,colorize=T,main=gsub("\\.","-",gsub("_.*="," => ",drug)))
  lines(x=c(0,1),y=c(0,1),lty=2)
  text(0.03,1,paste("AUC",try(round(unlist(performance(pred,"auc")@y.values),digits = 2),silent=T),sep="="))
  dev.off()
  
  result$Drug[counter]<-drug2
  result$Drug_Status[counter]<-drug.status
  result$Predictor_Genes[counter]<-paste(predictor.genes.0,collapse=",")
  result$Predictor_Genes_No[counter]<-length(predictor.genes.0)
  try(result$AUC[counter]<-unlist(performance(pred,"auc")@y.values),silent=T)
  try(result$Sensitivity[counter]<-unlist(performance(pred,"sens")@y.values)[which.max(unlist(performance(pred,"mat")@y.values))],silent=T)
  try(result$Specificity[counter]<-unlist(performance(pred,"spec")@y.values)[which.max(unlist(performance(pred,"mat")@y.values))],silent=T)
  try(result$PPV[counter]<-unlist(performance(pred,"ppv")@y.values)[which.max(unlist(performance(pred,"mat")@y.values))],silent=T)
  try(result$NPV[counter]<-unlist(performance(pred,"npv")@y.values)[which.max(unlist(performance(pred,"mat")@y.values))],silent=T)
  try(result$ACC[counter]<-unlist(performance(pred,"acc")@y.values)[which.max(unlist(performance(pred,"mat")@y.values))],silent=T)
  try(result$FPR[counter]<-unlist(performance(pred,"fpr")@y.values)[which.max(unlist(performance(pred,"mat")@y.values))],silent=T)
  counter=counter+1
  predictor.genes<-NA
  p.comb<-NA
  p.final<-NA
}
write.table(result,file="~/Projects/MBA/data/Prediction_h2o_rf_quantile_support_filter_GE_NUM/Predictors_Report.tsv",sep='\t',row.names = F)
h2o.shutdown(prompt = TRUE)
