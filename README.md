# NegBinSig-Test

This is a simple script that performs negative binomial and zero-inflated negative binomial regression.

The script is adjusted such that it can take in OTU table generated via QIIME from biom file as input and standard mapping/metadata file compatible with QIIME.

The script is run via command line using the Rscript command. To run the script, pass the command in following format:

Rscript nb_test_script.R both_vs_neither_filter.txt both_neither_mapfile.txt Both Neither Treatment ZINB_NB_Output_test_v2.txt 2

As seen from the command, the script takes in 7 commands. They are as follows:

1) OTU table generated via QIIME (which is called both_vs_neither_filter.txt in the above example)
2) QIIME compatible mapping file (which is called neither_mapfile.txt in the above example)
3) 

# where both_vs_neither_filter.txt = OTU table generated via QIIME
# neither_mapfile.txt = QIIME compatible mapping file
# Both = Category 1
# Neither = Category 2
# Treatment = Column name of the treatment group that needs to be compared. 
# ZINB_NB_Output_test_v2.txt = Name of output file
# 2 = No. of cores to use. More cores on machine, faster the analysis will complete
