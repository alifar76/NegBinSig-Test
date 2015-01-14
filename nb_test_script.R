start.time <- Sys.time()

library(pscl)
library(MASS)
library(foreach)
library(doMC)


zinb_nb_test <- function(both,trt){
    all_data <- foreach(i=2:(length(both)-1), .combine = rbind) %dopar% {
        final_vec <- c()
        result.zinb <- tryCatch(zeroinfl(both[ ,i] ~ trt | 1, data = both, dist = "negbin"),error=function(e) NULL)
        result.nb <- tryCatch(glm.nb(both[ ,i] ~ trt, data = both),error=function(e) NULL)
        if (!is.null(result.zinb) && !is.null(result.nb)){
            zinb.coeff <- exp(summary(result.zinb)$coefficients$count[2,1])
            nb.coeff <- exp(summary(result.nb)$coefficients[2,1])
            zinb.pval <- summary(result.zinb)$coefficients$count[2,4]
            nb.pval <- summary(result.nb)$coefficients[2,4]
            final_vec <- c(zinb.coeff,zinb.pval,nb.coeff,nb.pval)
        }
        if (is.null(result.zinb) && !is.null(result.nb)){
            nb.coeff <- exp(summary(result.nb)$coefficients[2,1])
            nb.pval <- summary(result.nb)$coefficients[2,4]
            final_vec <- c(NA,NA,nb.coeff,nb.pval)
        }
        if (!is.null(result.zinb) && is.null(result.nb)){
            zinb.coeff <- exp(summary(result.zinb)$coefficients$count[2,1])
            zinb.pval <- summary(result.zinb)$coefficients$count[2,4]
            final_vec <- c(zinb.coeff,zinb.pval,NA,NA)
        }
        if (is.null(result.zinb) && is.null(result.nb)){
            final_vec <- c(NA,NA,NA,NA)
        }
        shap_wilk_pval <- tryCatch(shapiro.test(both[,i])$p.value,error=function(e) NA)		# Significant p-value indicates data is not normally distributed.
        ttest <- t.test(both[,i] ~ trt, data=both)
		pval_ttest <- ttest$p.value
		estimate_tab <- ttest$estimate
		heading <- paste(gsub(" ","",strsplit(names(estimate_tab)[1],"mean in group")[[1]][2],fixed=TRUE),"_minus_",gsub(" ","",strsplit(names(estimate_tab)[2],"mean in group")[[1]][2],fixed=TRUE),"_mean",sep="")
		mean_diff <- estimate_tab[1][[1]] - estimate_tab[2][[1]]
      	final_vec <- c(final_vec,mean_diff,pval_ttest,shap_wilk_pval,heading)
        final_vec
    }
    return (all_data)
}


data_organizer_and_result <- function(otutable,mapfile,categ1,categ2,metavariable,outputname,num){
 
    registerDoMC(as.numeric(num))   #change the 4 to your number of CPU cores  
	MYdata <- read.table(otutable,header = T, sep = "\t", check.names = T, row.names =1, comment.char= "", skip =1)     # Ignore # at beginning of line. Skip first line which says "Converted from biom"
	MYmeta <- read.table(mapfile,header = T, sep = "\t", check.names = T, comment.char= "")
 	allcols <- length(colnames(MYdata))
    MYdata <- MYdata[,order(names(MYdata))]
    
    MYdata2 <- MYdata[,1:(allcols-1)]
    MYdata2 <- MYdata2[,order(names(MYdata2))]
    MYmeta <- MYmeta[order(MYmeta[,"X.SampleID"]), ]
    matrix1 <- MYdata2[c(as.factor(MYmeta[,"X.SampleID"][MYmeta[,metavariable]==categ1]))]
    matrix2 <- MYdata2[c(as.factor(MYmeta[,"X.SampleID"][MYmeta[,metavariable]==categ2]))]
    matrix3 <- cbind(matrix1,matrix2)
    mat.array <- t(matrix3)
    both <- merge(mat.array,MYmeta,by.x=0,by.y="X.SampleID")
    
    all_data <- zinb_nb_test(both,metavariable)
    zinb.qval <- p.adjust(all_data[,2], method = "fdr")
    nb.qval <-  p.adjust(all_data[,4], method = "fdr")
    ttest.qval <- p.adjust(all_data[,6], method = "fdr") 
    taxonomy <- MYdata[allcols]
    otuids <- colnames(both)[2:(length(colnames(both))-1)]
    taxlabels <- as.vector(taxonomy[,1])
    difflabel <- unique(all_data[,8])
    all_data <- cbind(otuids,all_data[,1:2],zinb.qval,all_data[,3:4],nb.qval,all_data[,5:6],ttest.qval,all_data[,7],taxlabels)
    colnames(all_data) <- c("OTU_IDs","ZINB_Coeff","ZINB_pval","ZINB_qval","NB_Coeff","NB_pval","NB_qval",difflabel,"ttest_pval","ttest_qval","Shapiro_Wilk_Normality_pvalue","taxonomy")
    write.table(as.matrix(all_data),file=outputname,sep="\t",append = TRUE,col.names=TRUE,row.names=FALSE,quote=FALSE)
}

argv <- commandArgs(TRUE)




data_organizer_and_result(argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7])

print (Sys.time() - start.time)
