This is a MATLAB implementation of causal computational models for gene regulatory network discovery (2015-16).

Undergraduate Students: Sahil Loomba and Parul Jain (both from CS departments at IIT Delhi)
Supervisor: Dr. Sumeet Agarwal, IIT Delhi

Instructions to the investigator:
1. Save your normalised gene expression dataset in a .dat file to ./data, in the form 'data_NAME.dat' where NAME is the (unique) name of your dataset. Ensure that this file consists of two variables, namely 'grn' and 'exp'.
2. To run pairwise tests, use master scripts in ./pairwise.
3. To run intrinsic graph estimation analysis, use master scripts in ./global.
4. View results in ./dat and plots in ./plots.

Note that a lot of the scripts right now have been suited for DREAM4 datasets. Thus, carefully go through the scripts to understand the changes you'll need to make for running it on other datasets.