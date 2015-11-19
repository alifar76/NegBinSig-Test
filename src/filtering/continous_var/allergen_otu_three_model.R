library(data.table)
library(pscl)
library(MASS)
library(foreach)
library(doMC)


allergen_otu_three_model <- function(otutable,mapfile,trtgrp,outputname){
	MYdata <- fread(otutable,header = T, sep = "\t", check.names = F, skip =1)
	MYmeta <- fread(mapfile,header=T,sep="\t",check.names=T)
	setcolorder(MYdata, sort(colnames(MYdata)))
	rowmeta <- sort(as.vector(MYmeta[,"#SampleID",with=FALSE])[[1]])
	MYmeta <- MYmeta[match(rowmeta, as.vector(MYmeta[,"#SampleID",with=FALSE])[[1]]),]
	all_data <- foreach(i=1:nrow(MYdata), .combine = rbind) %dopar% {
		otudat <- as.numeric(unlist(unname(MYdata[i,]))[2:(ncol(MYdata)-1)])
		trtvals <- as.vector(MYmeta[,trtgrp,with=FALSE])[[1]]
		result.pois <- tryCatch(glm(otudat ~ trtvals, family="poisson"),error=function(e) NA)
		result.nb <- tryCatch(glm.nb(otudat ~ trtvals),error=function(e) NA)
		result.zinb <- tryCatch(zeroinfl(otudat ~ trtvals | 1, dist = "negbin"),error=function(e) NA)
		out_pois <- tryCatch(summary(result.pois)$coefficients[2,],error=function(e) NA)
		out_nb <- tryCatch(summary(result.nb)$coefficients[2,],error=function(e) NA)
		out_zinb <- tryCatch(summary(result.zinb)$coefficients$count[2,],error=function(e) NA)
		if (is.na(out_pois)){
			out_pois <- rep(NA,4)
		} else {
			out_pois <- unname(out_pois)}
		if (is.na(out_nb)){
			out_nb <- rep(NA,4)
		} else {
			out_nb <- unname(out_nb)}
		if (is.na(out_zinb)){
			out_zinb <- rep(NA,4)
		} else {
			out_zinb <- unname(out_zinb)}
		bic.pois <- tryCatch(AIC(result.pois,k=log(length(otudat))),error=function(e) NA)											# Column 10	
		bic.nb <- tryCatch(AIC(result.nb,k=log(length(otudat))),error=function(e) NA)											# Column 11
		bic.zinb <- tryCatch(AIC(result.zinb,k=log(length(otudat))),error=function(e) NA)											# Column 12
		results <- c(MYdata[i,][[1]],out_pois,out_nb,out_zinb,bic.pois,bic.nb,bic.zinb,MYdata[i,][[ncol(MYdata)]])
		results
	}
	pois.qval <- p.adjust(all_data[,5],"fdr")
	nb.qval <- p.adjust(all_data[,9],"fdr")
	zinb.qval <- p.adjust(all_data[,13],"fdr")
	all_data <- cbind(all_data,pois.qval,nb.qval,zinb.qval)
	colnames(all_data) <- c("OTU_ID","Estimate (Poisson)","SE (Poisson)","z value (Poisson)","p-value (Poisson)","Estimate (NB)","SE (NB)","z value (NB)","p-value (NB)","Estimate (ZINB)","SE (ZINB)","z value (ZINB)","p-value (ZINB)","BIC (Poisson)","BIC (NB)", "BIC (ZINB)","taxonomy","Poisson (q-value)","NB (q-value)","ZINB (q-value)")
	write.table(as.matrix(all_data),file=outputname,sep="\t",append = TRUE,col.names=TRUE,row.names=FALSE,quote=FALSE)
}




otutable <- "ureca_sensitivity_248_samples.txt"
mapfile <- "allergen_1yr_musm_blag.txt"
nproc <- 4
trtgrp <- "musm_b_y1"
outputname <- "allergen_otu_impact_musm.txt"

trtgrp2 <- "blag_b_y1"
outputname2 <- "allergen_otu_impact_blag.txt"




registerDoMC(nproc)

sessionInfo()


allergen_otu_three_model(otutable,mapfile,trtgrp,outputname)
allergen_otu_three_model(otutable,mapfile,trtgrp2,outputname2)


