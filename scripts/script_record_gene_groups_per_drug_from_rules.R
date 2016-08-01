###################################
########Rules4#####################
###################################
load("~/Projects/MBA/data/rules4_proc.RData")

drugs<-unique(out4$rhs)
result<-data.frame(Drug=rep(NA,length(drugs)),
                   Drug_Status=rep(NA,length(drugs)),
                   Predictor_Genes=rep(NA,length(drugs)),
                   Predictor_Genes_No=rep(NA,length(drugs))
)
counter=1
for (drug in drugs){
  
  temp.rules<-out4[out4$rhs==drug,]
  ####Selection Criteria
  temp.rules<-temp.rules[with(temp.rules,order(-lift)),]#Sort by:
  temp.rules2<-temp.rules[1:200,]#All Genes
  ####
  
  temp.rules2$lhs<-gsub("=under|=over|=GAIN|=LOSS|=Mut","",temp.rules2$lhs)
  predictor.genes.0<-unique(temp.rules2$lhs)
  drug2<-unlist(strsplit(drug,"="))[1]
  drug.status<-unlist(strsplit(drug,"="))[2]
  result$Drug[counter]<-drug2
  result$Drug_Status[counter]<-drug.status
  result$Predictor_Genes[counter]<-paste(predictor.genes.0,collapse=",")
  result$Predictor_Genes_No[counter]<-length(predictor.genes.0)
  counter=counter+1
}
write.table(result,file="~/Projects/MBA/data/Genes_from_rules4_per_drug_top200lift.tsv",row.names = F,sep='\t')
###################################
########End of Rules4##############
###################################



