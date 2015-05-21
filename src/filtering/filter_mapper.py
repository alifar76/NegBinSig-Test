import os
from decimal import *

def zero_check(checkv,dictn):
	val = [Decimal(checkv)]
	result = [x/abs(x) for x in val if x !=0]
	if len(result) != 0:
		return dictn[str(result[0])]


def filter_results(iname,outname,qval,filtinds,qvals):
	infile = open(iname,'rU')
	outfile = open(outname,"w")
	grps = {}
	for line in infile:
		spline = line.strip().split("\t")
		if line.startswith('OTU_IDs'):
			# Model name, OTU ID, Mean Difference, Significant Group, Lowest BIC value, q-value for selected model, taxonomy
			outfile.write("Selected Model"+"\t"+spline[0]+"\t"+spline[12]+"\t"+"Significant Group"+"\t"+"BIC value of Selected Model (Lowest)"+"\t"+"q-value for selected model"+"\t"+spline[28]+"\t"+'\t'.join(spline[18:24])+"\n")		
			diffg = spline[12].split("_minus_")
			grps['1'] = diffg[0]																	# Treatment group 1 name
			grps['-1'] = diffg[1].replace("_mean","")												# Treatment group 2 name
		else:
			if spline[filtinds[0]:filtinds[1]] != ['NA','NA','NA']:									# Indices of BIC values of Poisson, NB and ZINB with outliers filtered/not filtered
				newc = filter(lambda a: a != 'NA', spline[filtinds[0]:filtinds[1]])					# Remove NAs from BIC values using lambda expression
				minbic = str(min([Decimal(e) for e in newc]))										# Convert the values to decimals and obtain the minimum value
				indexic = spline[filtinds[0]:filtinds[1]].index(minbic)								# Index of lowest value of BIC of filtered values
				if (indexic == 0) and (Decimal(spline[qvals[0]]) < qval):									# Filter OTUs with Poisson q-value less than 0.05
					outfile.write("Poisson"+"\t"+spline[0]+"\t"+spline[12]+"\t"+zero_check(spline[12],grps)+"\t"+spline[filtinds[0]]+"\t"+spline[qvals[0]]+"\t"+spline[28]+"\t"+'\t'.join(spline[18:24])+"\n")
				if (indexic == 1) and (Decimal(spline[qvals[1]]) < qval):									# Filter OTUs with NB q-value less than 0.05
					outfile.write("NB"+"\t"+spline[0]+"\t"+spline[12]+"\t"+zero_check(spline[12],grps)+"\t"+spline[filtinds[0]+1]+"\t"+spline[qvals[1]]+"\t"+spline[28]+"\t"+'\t'.join(spline[18:24])+"\n")
				if (indexic == 2) and (Decimal(spline[qvals[2]]) < qval):									# Filter OTUs with ZINB q-value less than 0.05
					outfile.write("ZINB"+"\t"+spline[0]+"\t"+spline[12]+"\t"+zero_check(spline[12],grps)+"\t"+spline[filtinds[0]+2]+"\t"+spline[qvals[2]]+"\t"+spline[28]+"\t"+'\t'.join(spline[18:24])+"\n")
	outfile.close()
	return outfile



# Infile name, output file name, q-value, indices of BIC values, indices of q-values

filtinds = [44,47]			# Indices of BIC values with outliers filtered
qvalfilt = [30,32,34] 		# Indices of q-values of 3 models with outliers filtered
filter_results("../ZINB_NB_Output_result.txt","results_filtered.txt",0.05,filtinds,qvalfilt)



notfilt = [38,41]			# Indices of BIC values with outliers not filtered
qvalnonfilt = [3,6,9]		# Indices of q-values of 3 models with outliers not filtered
filter_results("../ZINB_NB_Output_result.txt","results_not_filtered.txt",0.05,notfilt,qvalnonfilt)
