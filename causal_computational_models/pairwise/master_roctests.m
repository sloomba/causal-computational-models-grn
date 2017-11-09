%THIS SCRIPT WAS DEVELOPED IN THE EARLIER STAGES OF THE PROJECT. CAN BE IGNORED.
name = 'new'; %modify this to NAME of dataset
dataset = strcat('../data/data_',name); %set appropriate dataset name here
mkdir(strcat('../plots/plot_',name));
figname = strcat('../plots/plot_',name,'/'); %set appropriate destn plot folder here. comment if you wish to not save plots
resname = strcat('../data/result_',name,'_pw'); %set appropriate destn result mat name here
addpath('techniques');
addpath('../helpers');
load(dataset);
sizes = [100]; %[10,20,50,100];
times = [1000]; %[500,1000];
%bins = [2]; %[2,5,10,20]; %add Inf? uncomment if old-style quantised data is being used
%pw_out = {};
for i=1:length(sizes)
    siz = sizes(i);
    for j=1:length(times)
        tim = times(j);
        %for k=1:length(bins)
            %bin = bins(k);
            siz, tim%, bin
            fig = strcat(figname, 'GRN',num2str(siz),'nodes',num2str(tim),'timesteps');%,num2str(bin),'bins'); %comment if you wish to not save plots
            res = strcat(resname, 'GRN',num2str(siz),'nodes',num2str(tim),'timesteps');%,num2str(bin),'bins'); %comment if you wish to not save plots
            %if bin==Inf %uncomment this if-block if old-style quantised data is being used
            %    pw_out{siz,tim} = run_roctests(expn{siz,tim},grns{siz,tim},fig);
            %else
            %    pw_out{siz,tim,bin} = run_roctests(expn{siz,tim,bin},grns{siz,tim},fig);
            %end
		    pw_out = {};
            pw_out{siz,tim} = run_roctests(exp{siz,tim},grn{siz,tim},fig, res); %comment 'fig' if you wish to not save plots to an apt folder
%             vals = pw_out{siz,tim}; %uncomment this block if you want to do and save Correlation-CCM analysis
%             values = vals{6};
%             val_co = values{2};
%             val_cm = values{6};
%             plot_comp_co = [];
%             plot_comp_cm = [];
%             for l=1:size(val_co,1)
%                 for m=1:size(val_cm,1)
%                     if val_cm(m,2:3)==val_co(l,2:3)
%                         plot_comp_co = [plot_comp_co;val_co(l,end)];
%                         plot_comp_cm = [plot_comp_cm;val_cm(m,end)];
%                         break;
%                     end
%                 end
%             end
%             plot(plot_comp_co, plot_comp_cm,'o');
%             xlabel('Correlation of Link');
%             ylabel('CCM of Link');
%             title('Comparing Correlation and CCM');
%             for l=1:size(plot_comp_co,1)
%                 text(plot_comp_co(l), plot_comp_cm(l), strcat('(',num2str(val_co(l,1)),',',num2str(val_co(l,2)),')'));
%             end
%             print(strcat(fig, 'CompareCoCm'),'-dpng');
        %end
    end
end
%save(resname,'pw_out');       
