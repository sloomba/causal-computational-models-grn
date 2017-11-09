%Master script, to be run by investigator.
name = 'new'; %modify this to NAME of dataset
dataori = strcat('../data/data_',name);
dataset = strcat('../data/result_',name,'_pw'); %set appropriate dataset ("values") name here
mkdir(strcat('../plots/plot_',name));
figname = strcat('../plots/plot_',name,'/'); %set appropriate destn plot folder here. comment if you wish to not save plots
resname = strcat('../data/result_',name','_ige'); %set appropriate destn result mat name here
addpath('../helpers');
load(dataset,'pw_out'); %must have the variable 'out'
load(dataori);
sizes = [50, 100]; %[10,20,50,100];
times = [500,1000]; %[500,1000];
%bins = [2]; %[2,5,10,20]; %add Inf? uncomment if old-style quantised data is being used
method = containers.Map(1,'Random');
method(2) = 'Correlation';
method(3) = 'GrangerCausality';
method(4) = 'MutualInformation';
method(5) = 'TransferEntropy';
method(6) = 'ConvergentCrossMap';
ige_out = {};
for i=1:length(sizes)
    siz = sizes(i);
    for j=1:length(times)
        tim = times(j);
        %for k=1:length(bins)
            %bin = bins(k);
            siz, tim%, bin
            fig = strcat(figname, 'GRN',num2str(siz),'nodes',num2str(tim),'timesteps');%,num2str(bin),'bins'); %comment if you wish to not save plots
            vals = pw_out{siz,tim};
            vals = vals{7};
            %if bin==Inf %uncomment this if-block if old-style quantised data is being used
            %    ige_out{siz,tim} = ige_analysis(expn{siz,tim},grns{siz,tim},fig);
            %else
            %    ige_out{siz,tim,bin} = ige_analysis(expn{siz,tim,bin},grns{siz,tim},fig);
            %end
            for k=2:6
                metricmat = metriclist_to_metricmat(vals{k});
                ige_out{siz,tim,k} = ige_analysis(exp{siz,tim},grn{siz,tim},metricmat,strcat(fig,method(k))); %comment 'fig' if you wish to not save plots to an apt folder
            end
        %end
    end
end
save(resname,'ige_out');