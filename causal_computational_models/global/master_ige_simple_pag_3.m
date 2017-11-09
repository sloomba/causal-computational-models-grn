%Master script, to be run by investigator, for simple one-network datasets.
name = 'dream_10_3_correct'; %modify this to NAME of dataset
dataori = strcat('../data/data_',name);
dataset = strcat('../data/result_',name,'_pw'); %set appropriate dataset ("values") name here
mkdir(strcat('../plots/plot_',name));
figname = strcat('../plots/plot_',name,'/'); %set appropriate destn plot folder here. comment if you wish to not save plots
resname = strcat('../data/result_',name,'_pag'); %set appropriate destn result mat name here
addpath('../helpers');
load(dataset,'pw_out'); %must have the variable 'out'
load(dataori);
method = containers.Map(1,'Random');
method(2) = 'Correlation';
method(3) = 'GrangerCausality';
method(4) = 'MutualInformation';
method(5) = 'TransferEntropy';
method(6) = 'ConvergentCrossMap';
ige_out = {};
vals = pw_out{7};
for k=2:6
    metricmat = metriclist_to_metricmat(vals{k});
    ige_out{k} = ige_analysis_pag(grn,exp,metricmat,strcat(figname,method(k))); %comment 'fig' if you wish to not save plots to an apt folder
end
save(resname,'ige_out');