# NegBinSig-Test

Background
------

This is a simple script that performs negative binomial and zero-inflated negative binomial regression.

The script is adjusted such that it can take in OTU table generated via QIIME from biom file as input and standard mapping/metadata file compatible with QIIME.

Presently, the script can only perform a single category comparison for variables. For example, if the metadata have two
variables such as diet and antibiotic exposure, the script will have to be run seperately for each variable. A joint model
including both explanatory variables (i.e., diet and antibiotic exposure) cannot be currently calculated.

Furthermore, the script can only work with variables that that have two levels. For example, if the 
explanatory variable is temperature, it must only contain two levels such as High and Low for the script to run. 

Work is in progress to include comparisons for more than 2 levels and to create models in which multiple combinations of 
variables can be included.

Running the script
------

The script is run via command line using the Rscript command. To run the script, pass the command in following format:

```Rscript nb_test_script.R both_vs_neither_filter.txt both_neither_mapfile.txt Both Neither Treatment ZINB_NB_Output_test_v2.txt 2```

As seen from the command, the script takes in 7 commands. They are as follows:

1) OTU table generated via QIIME (which is called **both_vs_neither_filter.txt** in the above example)

2) QIIME compatible mapping file (which is called **neither_mapfile.txt** in the above example)

3) Level 1 of the category being compared (which is called **Both** in the above example)

4) Level 2 of the category being compared (which is called **Neither** in the above example)

5) Column name of the category being compared as labelled in the mapping file (which is called **Treatment** in the above example)

6) Output file contaiing result (which is called **ZINB_NB_Output_test_v2.txt** in the above example)

7) No. of cores to use. More cores on machine, faster the analysis will complete (which is **2** in the above example)
