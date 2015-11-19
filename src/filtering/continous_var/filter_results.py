import os
from decimal import *



# Pois, NB, ZINB
bicvals = [13,16]		# 13:15 is range of bic values
qvals = [17,18,19]		
estimate = [1,5,9]		

# Indices of more variables
taxon = 16
otu_id = 0
qval = 0.05


def file_filterer(infile,ofile,taxon,otu_id,qval,bicvals,qvals,estimate):
	outfile = open(ofile,"w")
	infile = open(infile,'rU')
	for line in infile:
		if line.startswith("OTU_ID"):
			outfile.write("Model"+"\t"+"OTU ID"+"\t"+"Estimate"+"\t"+"q-value"+"\t"+"taxonomy"+"\n")
		else:
			spline = line.strip().split("\t")
			if spline[bicvals[0]:bicvals[1]] != ['NA','NA','NA']:
				newc = filter(lambda a: a != 'NA', spline[bicvals[0]:bicvals[1]])					# Remove NAs from BIC values using lambda expression
				minbic = str(min([Decimal(e) for e in newc]))										# Convert the values to decimals and obtain the minimum value
				indexic = spline[bicvals[0]:bicvals[1]].index(minbic)								# Index of lowest value of BIC of filtered values
				if spline[qvals[2]] != 'NaN':
					if (indexic == 0) and (Decimal(spline[qvals[0]]) < qval):
						outfile.write("Poisson"+"\t"+spline[otu_id]+"\t"+spline[estimate[0]]+"\t"+spline[qvals[0]]+"\t"+spline[taxon]+"\n")
					if (indexic == 1) and (Decimal(spline[qvals[1]]) < qval):
						outfile.write("NB"+"\t"+spline[otu_id]+"\t"+spline[estimate[1]]+"\t"+spline[qvals[1]]+"\t"+spline[taxon]+"\n")
					if (indexic == 2) and (Decimal(spline[qvals[2]]) < qval):
						outfile.write("ZINB"+"\t"+spline[otu_id]+"\t"+spline[estimate[2]]+"\t"+spline[qvals[2]]+"\t"+spline[taxon]+"\n")
	outfile.close()
	return outfile


file_filterer("allergen_otu_impact_musm.txt","musm_filtered_otus.txt",taxon,otu_id,qval,bicvals,qvals,estimate)
file_filterer("allergen_otu_impact_blag.txt","blag_filtered_otus.txt",taxon,otu_id,qval,bicvals,qvals,estimate)






