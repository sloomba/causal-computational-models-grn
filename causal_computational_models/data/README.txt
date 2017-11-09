This folder contains .dat files.
Those labeled as 'data_NAME.mat' are the GRN datasets, consisting of two variables, namely 'grn' and 'exp'.
These datasets could either be simple (consisting of only one sample network: such as 'data_dream.mat') or composite (consisting of multiple sample networks in a cell structure, indexed by {size,timesteps}: such as 'data_new.mat'). Note: 'data_old.mat' must not be used, not a good dataset. 
Those labeled as 'result_NAME_pw.mat' are the parwise analysis results, consisting of one (cell) variable 'pw_out' = {precision, recall, fscore, fpr, tpr, values, unsorted_values}.
Those labeled as 'result_NAME_ige.mat' are the intrinsic graph estimation analysis results, consisting of one (cell) variable 'ige_out' = {p, r, f, a, adj_found, roc, auroc, precision, recall, fscore}.
Note that you need only worry about DREAM4 data files.