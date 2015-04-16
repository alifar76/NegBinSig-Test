# NegBinSig-Test

Background
------

This is a simple script that performs negative binomial and zero-inflated negative binomial regression.

The script is adjusted such that it can take in any OTU table file, generated via [QIIME 1.8.0 (stable public release)](http://qiime.org/), (in tab-delimited format) as input and standard mapping/metadata file compatible with QIIME.

Presently, the script can only perform a single category comparison for variables. For example, if the metadata have two
variables such as diet and antibiotic exposure, the script will have to be run seperately for each variable. A joint model
including both explanatory variables (i.e., diet and antibiotic exposure) cannot be currently calculated.

Furthermore, the script can only work with variables that that have two levels. For example, if the 
explanatory variable is temperature, it must only contain two levels such as High and Low for the script to run. 

Work is in progress to include comparisons for more than 2 levels and to create models in which multiple combinations of explanatory variables can be included.

Running the script
------

There are 3 scripts in the folder src. The one to use is called ```nb_regression_outlier_filtering.R```. The script is run via command line using the Rscript command (in terminal). To run the script, pass the command in following format:

```Rscript nb_regression_outlier_filtering.R high_vs_low_otu_table.txt high_low_mapfile.txt High Low Treatment ZINB_NB_Output_result.txt 2```

As seen from the command, the script takes in 7 commands. They are as follows:

1) OTU table generated via QIIME (which is called **high_vs_low_otu_table.txt** in the above example)

2) QIIME compatible mapping file (which is called **high_low_mapfile.txt** in the above example)

3) Level 1 of the category being compared (which is called **High** in the above example)

4) Level 2 of the category being compared (which is called **Low** in the above example)

5) Column name of the category being compared as labelled in the mapping file (which is called **Treatment** in the above example)

6) Output file contaiing result (which is called **ZINB_NB_Output_result.txt** in the above example)

7) No. of cores to use. More cores on machine, faster the analysis will complete (which is **2** in the above example)

Please ensure that all the 8 arguments are provided, in the correct order and format. Otherwise, the script will crash and cause problems.

Input file format
------

Input of file format should be one compatabile with QIIME. However, please ensure that the sample IDs are not numeric. That is, the sample IDs should not be like: 1560.1, 1561.1, 1559.1, etc. If such is the case, please slightly modify the sample IDs in both the mapping file and OTU table by adding any alphabet. So, for example, sample ID 1560.1 will become p1560.1.

Also, please make sure that the mapping file has the same number of samples as the OTU tables, having the same sample IDs. If mapping file has more or less sample IDs than the samples in the OTU table, the script will crash.

Output Explained
------

The output of the script contains information for all of the OTUs tested. Currently, there are 51 columns in the output file (called **ZINB_NB_Output_result.txt** in the example above) as generated via this script. The columns and their descriptions of the output file are as follows:

1) **OTU_ID**: Indicates the OTU ID

2) **Poiss_Coeff**: Indicates the expoentiated regression coefficient estimated when fitting the Poisson model on the OTU count data. This value is the is the multiplicative change in OTU abundance, comparing one level of a specific treatment group to another. It is to be interpreted in a similar way to the NB_Coeff value, which is explained further.

3) **Poiss_pval**: p-value of estimated Poisson Coeff

4) **Poiss_qval**: q-value of the estimated Poisson Coeff

5) **NB_Coeff**: Indicates the exponentiated regression coefficient for the regular negative binomial model. To elaborate more,
in our example dataset, OTU_26 has a NB_Coeff value of 1.329890281. This means that the abundance of OTU_26 is 1.329890281 times higher in the Low group of treatment compared to the High group of treatment. To find out which group is the base group, we look at the column called **High_minus_Low_mean** for our example dataset. (The name of this group will change from mapping file to mapping file as explained later in detail) In our example dataset, the **High_minus_Low_mean** column has the value of -176.8338789. This suggests that the mean OTU_26 is higher in Low group when compared to High group.

6) **NB_pval**: p-value of the estimated NB_Coeff

7) **NB_qval**: q-value of the estimated NB_Coeff

8) **ZINB_Coeff**: Indicates the exponentiated regression coefficient for the zero inflated negative binomial model. This value is the multiplicative change in OTU abundance, comparing one level of a specific treatment group to another. It is to be interpreted in a similar way as NB_Coeff value, which is explained earlier. For many OTUs, this value may turn out to be NA and it's okay if it does. It simply reflects that ZINB is not a good model to fit to data (mostly due to convergence issues).

9) **ZINB_pval**: p-value of the estimated ZINB Coeff

10) **ZINB qval**: q-value of the estimated ZINB Coeff

11) **mean in group High**: The mean of a specific OTU in a samples belonging to specific one treatment level. (The name of this column will change base on user defined mapping file).

12) **mean in group Low**: The mean of a specific OTU in a samples belonging to specific the other treatment level. (The name of this column will change base on user defined mapping file).

13) **High_minus_Low_mean**: The name of this column is specific to each dataset. In our metadata file, we had a column called Treatment which had two levels: High and Low. The ordering of the names of the treatment levels in this column will vary from mapping file to mapping file. However, they will always be consistent such that if the value is positive, then the first level as suggested by column name has higher mean than second one (High in our example). And if the value is negative, then the second level as suggested by column name has a higher mean than first one (Low in our example).

14) **ttest_pval**: p-value of the t-test.

15) **ttest_qval**: q-value of the t-test.

16) **KW_pval**: Kruskal-Wallis pvalue.

17) **KW_qval**: Kruskal-Wallis qvalue.

18) **NB_Coeff_Estimate_Error**: Indicates (with yes/no) whether there was in error in convergence in estimating the negative binomial regression coefficient.

19) **# of 0's in High**: Indicates the total # of zeroes in one treatment group (such as High).

20) **# of 0's in Low**: Indicates the total # of zeroes in the other treatment group (such as Low).

21) **# of non-zeroes in High**: Indicates the total number of non-zero entries in one treatment group (such as High).

22) **# of non-zeroes in Low**: Indicates the total number of non-zero entries in other treatment group (such as Low).

23) **Total count in High**: Sum of all values in one treatment group (such as High).

24) **Total count in Low**: Sum of all values in other treatment group (such as Low).

25) **mean_otu**: Mean across all treatment groups (both High and Low).

26) **variance_otu**: Variance across all treatment groups (both High and Low).

27) **var/mean ratio**: Variance to mean ratio across all treatment groups (such as High and Low).

28) **Shapiro_Wilk_Normality_pvalue**: Indicates whether the data is normally distributed or not, informing us about the validity of using the t-test. A significant p-value in this column indicates that data are not normally distributed and t-test may not be that appropriate. 

29) **taxonomy**: Indicates the taxonomy/lineage of the specific OTU.

30) **pois_filt_pval**: p-value of the Poisson model with outlier(s) filtered. An outlier for a specific OTU is defined as that count value for the OTU which is greater than 5 times the inter-quartile range of the OTU across all samples.

31) **pois_filt_qval**: q-value of Poisson model with outlier(s) filtered.

32) **nb_filt_pval**: p-value of negative binomial model with outlier(s) filtered. Outlier defined as above.

33) **nb_filt_qval**: q-value of negative binomial model with outlier(s) filtered.

34) **zinb_filt_pval**: p-value of zero inflated negative binomial model with outlier(s) filtered. Outlier defined as above.

35) **zinb_filt_qval**: q-value of zero inflated negative binomial model with outlier(s) filtered.

36) **aic.pois**: [Akaike information criterion (AIC)] (http://en.wikipedia.org/wiki/Akaike_information_criterion) value for Poisson model.

37) **aic.nb**: [Akaike information criterion (AIC)] (http://en.wikipedia.org/wiki/Akaike_information_criterion) value for negative binomial model.

38) **aic.zinb**: [Akaike information criterion (AIC)] (http://en.wikipedia.org/wiki/Akaike_information_criterion) value for zero inflated negative binomial model.

39) **bic.pois**: [Bayesian information criterion (BIC)] (http://en.wikipedia.org/wiki/Bayesian_information_criterion) value for Poisson model.

40) **bic.nb**: [Bayesian information criterion (BIC)] (http://en.wikipedia.org/wiki/Bayesian_information_criterion) value for negative binomial model.

41) **bic.zinb**: [Bayesian information criterion (BIC)] (http://en.wikipedia.org/wiki/Bayesian_information_criterion) value for zero inflated negative binomial model.

42) **aic.filt.pois**: [Akaike information criterion (AIC)] (http://en.wikipedia.org/wiki/Akaike_information_criterion) value for Poisson model with outlier(s) filtered.

43) **aic.filt.nb**: [Akaike information criterion (AIC)] (http://en.wikipedia.org/wiki/Akaike_information_criterion) value for negative binomial model with outlier(s) filtered.

44) **aic.filt.zinb**: [Akaike information criterion (AIC)] (http://en.wikipedia.org/wiki/Akaike_information_criterion) value for zero inflated negative binomial model with outlier(s) filtered.

45) **bic.filt.pois**: [Bayesian information criterion (BIC)] (http://en.wikipedia.org/wiki/Bayesian_information_criterion) value for Poisson model with outlier(s) filtered.

46) **bic.filt.nb**: [Bayesian information criterion (BIC)] (http://en.wikipedia.org/wiki/Bayesian_information_criterion) value for negative binomial model with outlier(s) filtered.

47) **bic.filt.zinb**: [Bayesian information criterion (BIC)] (http://en.wikipedia.org/wiki/Bayesian_information_criterion) value for zero inflated negative binomial model with outlier(s) filtered.

48) **aic.nonfilt.best**: Column indicating which model is best (based on lowest AIC value) for a given OTU with no outliers filtered.

49) **bic.nonfilt.best**: Column indicating which model is best (based on lowest BIC value) for a given OTU with no outliers filtered.

50) **aic.filt.best**: Column indicating which model is best (based on lowest AIC value) for a given OTU with outliers filtered.

51) **bic.filt.best**: Column indicating which model is best (based on lowest BIC value) for a given OTU with outliers filtered.
