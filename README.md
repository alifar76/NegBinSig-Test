# NegBinSig-Test

Background
------

This is a simple script that performs negative binomial and zero-inflated negative binomial regression.

The script is adjusted such that it can take in any OTU table file, generated via QIIME, (in tab-delimited format) as input and standard mapping/metadata file compatible with QIIME.

Presently, the script can only perform a single category comparison for variables. For example, if the metadata have two
variables such as diet and antibiotic exposure, the script will have to be run seperately for each variable. A joint model
including both explanatory variables (i.e., diet and antibiotic exposure) cannot be currently calculated.

Furthermore, the script can only work with variables that that have two levels. For example, if the 
explanatory variable is temperature, it must only contain two levels such as High and Low for the script to run. 

Work is in progress to include comparisons for more than 2 levels and to create models in which multiple combinations of 
variables can be included.

Running the script
------

The script is run via command line using the Rscript command (in terminal). To run the script, pass the command in following format:

```Rscript nb_test_script.R both_vs_neither_filter.txt both_neither_mapfile.txt Both Neither Treatment ZINB_NB_Output_test_v2.txt 2```

As seen from the command, the script takes in 7 commands. They are as follows:

1) OTU table generated via QIIME (which is called **both_vs_neither_filter.txt** in the above example)

2) QIIME compatible mapping file (which is called **neither_mapfile.txt** in the above example)

3) Level 1 of the category being compared (which is called **Both** in the above example)

4) Level 2 of the category being compared (which is called **Neither** in the above example)

5) Column name of the category being compared as labelled in the mapping file (which is called **Treatment** in the above example)

6) Output file contaiing result (which is called **ZINB_NB_Output_test_v2.txt** in the above example)

7) No. of cores to use. More cores on machine, faster the analysis will complete (which is **2** in the above example)

Output Explained
------

The output of the script contains information for all of the OTUs tested. Currently, there are 12 columns in the output file. They are as follows:

1) **OTU_ID**: Indicates the OTU ID

2) **ZINB_Coeff**: Indicates the exponentiated regression coefficient for the zero inflated negative binomial model. This value is the multiplicative change in OTU abundance, comparing one level of a specific treatment group to another. It is to be interpreted in a similar way to the NB_Coeff value, which is explained further. For many OTUs, this value may turn out to be NA and it's okay if it does. It simply reflects that ZINB is not a good model to fit to data (mostly due to convergence issues).

3) **ZINB_pval**: p-value of the estimated ZINB Coeff

4) **ZINB qval**: q-value of the estimated ZINB Coeff

5) **NB_Coeff**: Indicates the exponentiated regression coefficient for the regular negative binomial model. To elaborate more,
in our example dataset, OTU_26 has a NB_Coeff value of 1.685675928. This means that the abundance of OTU_26 is 1.685675928 times higher in the Neither group of treatment compared to the Both group of treatment. To find out which group is the base group, we look at the column called **Both_minus_Neither_mean**. In our example dataset, the **Both_minus_Neither_mean** column has the value of -2324.903642. This suggests that the mean OTU_26 is higher in Neither group when compared to Both group.

6) **NB_pval**: p-value of the estimated NB_Coeff

7) **NB_qval**: q-value of the estimated NB_Coeff

8) **Both_minus_Neither_mean**: The name of this column is specific to each dataset. In our metadata file, we had a column called Treatment which had two levels: Both and Neither. The ordering of the names of the treatment levels in this column will vary from mapping file to mapping file. However, they will be always be consistent such that if the value is positive, then the first level as suggested by column name has higher mean than second one (Both in our example). And if the value is negative, then the second level as suggested by column name has a higher mean than first one (Neither in our example).

9) ttest_pval: p-value of the t-test

10) ttest_qval: q-value of the t-test

11) **Shapiro_Wilk_Normality_pvalue**: Indicates whether the data is normally distributed or not, informing us about the valdiity of using the t-test. A significant p-value in this column indicates that data are not normally distributed and t-test may not be that appropriate. 

12) **taxonomy**: Indicates the taxonomy/lineage of the specific OTU.
