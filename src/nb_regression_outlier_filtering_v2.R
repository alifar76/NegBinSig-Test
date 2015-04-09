#Rscript nb_regression_outlier_filtering_v2.R high_vs_low_otu_table.txt high_low_mapfile.txt High Low Treatment ZINB_NB_Output_result.txt 2

start.time <- Sys.time()

library(pscl)
library(MASS)
library(foreach)
library(doMC)

 
data_load <- function(MYdata,mapfile,treatment,trt1,trt2){    
  colnames(MYdata)
  MYmeta <- read.table(mapfile,header = T, sep = "\t", check.names = F, comment.char= "") #change Group header to Treatment
  colnames(MYmeta)
  allcols <- length(colnames(MYdata))
  MYdata <- MYdata[,order(names(MYdata))]
  MYdata2 <- MYdata[,1:(allcols-1)]
  MYdata2 <- MYdata2[,order(names(MYdata2))]
  MYmeta <- MYmeta[order(MYmeta[,"#SampleID"]), ]
  matrix1 <- MYdata2[c(as.factor(MYmeta[,"#SampleID"][MYmeta[,treatment]==trt1]))] # change whichever group you are testing
  matrix2 <- MYdata2[c(as.factor(MYmeta[,"#SampleID"][MYmeta[,treatment]==trt2]))] # change whichever group you are testing
  matrix3 <- cbind(matrix1,matrix2)
  mat.array <- t(matrix3)
  both <- merge(mat.array,MYmeta,by.x=0,by.y="#SampleID")
}
 

zinb_nb_test <- function(both,MYdata,trt,categ1,categ2){
  all_data <- foreach(i=2:(length(rownames(MYdata))+1), .combine = rbind) %dopar% {
    final_vec <- c()
    formula1 <- as.formula(paste("both[,i] ~ ",trt," | 1",sep=""))
    formula2 <- as.formula(paste("both[,i] ~ ",trt,sep=""))
    result.zinb <- tryCatch(zeroinfl(formula1, data = both, dist = "negbin"),error=function(e) NA)
    result.nb <- tryCatch(glm.nb(formula2, data = both),error=function(e) NA)
    zinb.coeff <- tryCatch(exp(summary(result.zinb)$coefficients$count[2,1]),error=function(e) NA)
    nb.coeff <- tryCatch(exp(summary(result.nb)$coefficients[2,1]),error=function(e) NA)
    zinb.pval <- tryCatch(summary(result.zinb)$coefficients$count[2,4],error=function(e) NA)
    nb.pval <- tryCatch(summary(result.nb)$coefficients[2,4],error=function(e) NA)
    final_vec <- c(zinb.coeff,zinb.pval,nb.coeff,nb.pval)
    shap_wilk_pval <- tryCatch(shapiro.test(both[,i])$p.value,error=function(e) NA)       # Significant p-value indicates data is not normally distributed.
    pval_ttest <- tryCatch(t.test(formula2, data=both)$p.value,error=function(e) NA)
    estimate_tab <- tryCatch(t.test(formula2, data=both)$estimate,error=function(e) NA)
    heading <- paste(gsub(" ","",strsplit(names(estimate_tab)[1],"mean in group")[[1]][2],fixed=TRUE),"_minus_",gsub(" ","",strsplit(names(estimate_tab)[2],"mean in group")[[1]][2],fixed=TRUE),"_mean",sep="")
    mean_diff <- tryCatch((estimate_tab[1][[1]] - estimate_tab[2][[1]]),error=function(e) NA)
    kwtest <- tryCatch(kruskal.test(formula2,data=both)$p.value,error=function(e) NA)
    warn.nb <- tryCatch(glm.nb(formula2, data = both),error=function(e) NA,warning=function(w) w)
    valwarn.nb <- ifelse(class(warn.nb)[1] == "simpleWarning", 'yes', 'no')
    trt1vals <- both[,i][which(as.vector(both[,trt]) == categ1)]
    trt2vals <- both[,i][which(as.vector(both[,trt]) == categ2)]
    zerotrt1 <- sum(trt1vals == 0)
    zerotrt2 <- sum(trt2vals == 0)
    nonzerotrt1 <- length(trt1vals) - zerotrt1
    nonzerotrt2 <- length(trt2vals) - zerotrt2
    totaltrt1 <- sum(trt1vals)
    totaltrt2 <- sum(trt2vals)
    mean.otu <- mean(both[,i])
    var.otu <- var(both[,i])
    var.mean.ratio <- var.otu/mean.otu
    
    newd <- both[,i][-(which(both[,i] > 5*IQR(both[,i])))]
	 treatment <- as.vector(both[,trt])
	 newmeta <- treatment[-(which(both[,i] > 5*IQR(both[,i])))]    # Select values greater than 5 times the IQR
	 newpval_nb <- tryCatch(summary(glm.nb(newd ~ newmeta))$coefficients[2,4],error=function(e) NA)

    final_vec <- c(final_vec,mean_diff,pval_ttest,shap_wilk_pval,heading,estimate_tab[1][[1]],names(estimate_tab)[1],estimate_tab[2][[1]],names(estimate_tab)[2],kwtest,valwarn.nb,zerotrt1,zerotrt2,nonzerotrt1,nonzerotrt2,totaltrt1,totaltrt2,mean.otu,var.otu,var.mean.ratio,newpval_nb)
    final_vec
  }
  return (all_data)
}
 

final_steps <- function(otutable,mapfile,categ1,categ2,trt,outputname){
  MYdata <- read.table(otutable,header = T, sep = "\t", check.names = F, row.names =1, comment.char= "", skip =1,quote="")
  both <- data_load(MYdata,mapfile,trt,categ1,categ2)
  all_data <- zinb_nb_test(both,MYdata,trt,categ1,categ2)
  allcols <- length(colnames(MYdata))
  zinb.qval <- p.adjust(all_data[,2], method = "fdr")
  nb.qval <-  p.adjust(all_data[,4], method = "fdr")
  ttest.qval <- p.adjust(all_data[,6], method = "fdr")
  kw.qval <- p.adjust(all_data[,13], method = "fdr")
  outlier_filt_qval <- p.adjust(all_data[,24], method = "fdr")
  taxonomy <- MYdata[allcols]
  otuids <- colnames(both)[2:(length(colnames(both))-1)]
  taxlabels <- as.vector(taxonomy[,1])
  difflabel <- unique(all_data[,8])
  mean1_head <- unique(all_data[,10])
  mean2_head <- unique(all_data[,12])
 zrtrt1 <- paste("# of 0's in ",categ1,sep="")
 zrtrt2 <- paste("# of 0's in ",categ2,sep="")
 nzrtrt1 <- paste("# of non-zeroes in ",categ1,sep="")
 nzrtrt2 <- paste("# of non-zeroes in ",categ2,sep="")
 tottrt1 <- paste("Total count in ",categ1,sep="")
 totttrt2 <- paste("Total count in ",categ2,sep="")
  all_data <- cbind(otuids,all_data[,1:2],zinb.qval,all_data[,3:4],nb.qval,all_data[,9],all_data[,11],all_data[,5:6],ttest.qval,all_data[,13],kw.qval,all_data[,14:23],all_data[,7],taxlabels,all_data[,24],outlier_filt_qval)
  colnames(all_data) <- c("OTU_IDs","ZINB_Coeff","ZINB_pval","ZINB_qval","NB_Coeff","NB_pval","NB_qval",mean1_head,mean2_head,difflabel,"ttest_pval","ttest_qval","KW_pval","KW_qval","NB_Coeff_Estimate_Error",zrtrt1,zrtrt2,nzrtrt1,nzrtrt2,tottrt1,totttrt2,"mean_otu","variance_otu","var/mean ratio","Shapiro_Wilk_Normality_pvalue","taxonomy","outlier_nbpval","outlier_nbqval") #change Difference to groups being tested
  write.table(as.matrix(all_data),file=outputname,sep="\t",append = TRUE,col.names=TRUE,row.names=FALSE,quote=FALSE)
}



argv <- commandArgs(TRUE)

registerDoMC(as.numeric(argv[7]))   #change the 4 to your number of CPU cores
final_steps(argv[1],argv[2],argv[3],argv[4],argv[5],argv[6])
print (Sys.time() - start.time)

