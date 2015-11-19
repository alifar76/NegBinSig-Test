import re

hierarchy = ["k__","p__","c__","o__","f__","g__","s__"]

def level_namer(num,hierarchy,test,listn):
	if test == 'Unassigned':
		return test
	else:
		sp = hierarchy[num]			# -1
		grp = re.search(r'(%s.*)'%sp,test)
		if grp:
			level = grp.group(1).split(sp)[1].split(":")[0]
			if level == '':
				num = num - 1
				return level_namer(num,hierarchy,test,listn)
			else:
				listn.append(level)
				if len(listn) == 2:
					first_name = listn[1].split(":")[0].replace('[','').replace(']','')
					second_name = listn[0].replace('[','').replace(']','')
					return first_name, second_name 
				else:
					num = num - 1
					return level_namer(num,hierarchy,test,listn)
		else:
			num = num - 1
			return level_namer(num,hierarchy,test,listn)

def main_function(f,outn):
	otudict = {}
	infile = open(f,'rU')
	for line in infile:
		if not line.startswith("Selected Model"):
			spline = line.strip().split("\t")
			otudict[spline[1]] = spline[6]
	outfile = open(outn,"w")
	for x in otudict:
		annotations = level_namer(-1,hierarchy,otudict[x],[])
		if annotations == 'Unassigned':
			outfile.write(x+"\t"+annotations+"\n")
		else:
			outfile.write(x+"\t"+'\t'.join(annotations)+"\n")
	outfile.close()
	return outfile

main_function("results_not_filtered.txt","OTU_annotations.txt")

#tests = ["k__Bacteria: p__Firmicutes: c__Clostridia: o__Clostridiales: f__: g__: s__",
#"k__Bacteria: p__Proteobacteria: c__Betaproteobacteria: o__Burkholderiales: f__Comamonadaceae: g__: s__",
#"k__Bacteria: p__Armatimonadetes: c__[Fimbriimonadia]: o__[Fimbriimonadales]: f__[Fimbriimonadaceae]: g__Fimbriimonas: s__"]
#for test in tests:
#	print level_namer(-1,hierarchy,test,[])
