%Master script, to be run by investigator, for simple one-network datasets.
name = 'dream_100_5_correct'; %modify this to NAME of dataset, 10/100?
dataset = strcat('../data/data_',name); %set appropriate dataset name here
mkdir(strcat('../plots/plot_',name));
figname = strcat('../plots/plot_',name,'/'); %set appropriate destn plot folder here. comment if you wish to not save plots
resname = strcat('../data/result_',name,'_pw'); %set appropriate destn result mat name here
addpath('techniques');
addpath('../helpers');
load(dataset);
pw_out = run_roctests(exp, grn, figname, resname); %comment 'fig' if you wish to not save plots to an apt folder
% save(resname,'pw_out');
            
% out = {precision, recall, fscore, fpr, tpr, values, unsorted_values};