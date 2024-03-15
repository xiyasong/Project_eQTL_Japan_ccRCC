
############################################### main KM analysis for gene-----------------------------

rm(list=ls())

library(survival)
setwd("/Users/xiyas/eQTL-continue/eQTL-continue-survival-analysis")
source('Cheng_toolbox_beta_gene.R')

cancerType<-"ccRCC"

#path_raw<-"E:\\L1000_2021\\kidney_drug_reposition\\coexp_network\\based_all_genes\\TPM_1\\Drug_reposition\\Cor_SHvsCP_HA1E\\KM_anaylysis\\target_gene\\"
#path_raw<-paste0(path_raw,"JAP\\")
path_raw<-"/Users/xiyas/eQTL-continue/eQTL-continue-survival-analysis/"
setwd(path_raw)
exp<-as.matrix(read.table(paste0(cancerType,"_tpm_target_gene_exp.txt"),header=T,sep="\t"))
exp
gene_list<-as.matrix(colnames(exp)[8:dim(exp)[2]])

dataDir=path_raw
path_raw
path_out_pdf<-paste0(path_raw,"KM_pdf/",sep = "")
path_out_pdf
dir.create(path_out_pdf)
setwd(path_out_pdf)

p_list<-NULL
cutoff_list<-NULL
coef_list<-NULL
for (j in 1:length(gene_list)){
  print(j)
  output<-Cheng_generateSurvInputfromTCGA(gene_list[j,1],cancerType,dataDir)
  setwd(path_out_pdf)
  result<-Cheng_generateKMplot(output,outFile=paste0(cancerType,"_",gene_list[j,1]))
  log_rank_p<-as.matrix(result$logRankP)
  cut_off<-as.matrix(result$EXPcut)
  coef<-as.matrix(result$coef)
  
  p_list<-rbind(p_list,log_rank_p)
  cutoff_list<-rbind(cutoff_list,cut_off)
  coef_list<-rbind(coef_list,coef)
  rm(log_rank_p,cut_off,coef)
}

fdr_list<-p.adjust(p_list,method="BH",n=length(p_list))
final_result<-cbind(as.matrix(gene_list),coef_list,cutoff_list,p_list,fdr_list)
colnames(final_result)<-c("gene","coef","cutoff","p","fdr")
write.table(final_result,file="coef_cutoff_list.txt",sep="\t",row.names=F,col.names=T,quote=F)

final_result <- as.data.frame(final_result)
KM_0.05 <- filter(final_result,p<= 0.05)
fdr_0.05 <- filter(final_result,fdr <=0.05)
KM_0.05_genes <- KM_0.05$gene
