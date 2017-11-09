%Master script, to be run by investigator, for simple one-network datasets.
name = 'new_GRN100nodes500timesteps'; %modify this to NAME of dataset
dataori = '../data/data_new';	%strcat('../data/data_',name);
dataset = strcat('../data/result_',name,'_pw'); %set appropriate dataset ("values") name here
mkdir(strcat('../plots/plot_',name));
figname = strcat('../plots/plot_',name,'/'); %set appropriate destn plot folder here. comment if you wish to not save plots
resname = strcat('../data/result_',name,'_ige'); %set appropriate destn result mat name here
addpath('../helpers');
load(dataset,'pw_out'); %must have the variable 'out'
load(dataori);
grn = grn{100, 500};
exp = exp{100, 500};
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
    ige_out{k} = ige_analysis(grn,exp,metricmat,strcat(figname,method(k))); %comment 'fig' if you wish to not save plots to an apt folder
end
save(resname,'ige_out');