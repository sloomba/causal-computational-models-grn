function [ out ] = topk_edge_analysis( num_nodes, suffix, lagtest )
%OUTPUT_AUROCS Outputs AUROCS for PW methods on DREAM4 datasets
%num_nodes = {10, 100}
%suffix = {'pw','ige','pag','jug'} %only pw right now
%lagtest = {'', 'CO', 'GC', 'CCM'}
%Example usage: output_aurocs(10,'pw','CCM');
if nargin<3
    lagtest = ''; %regular PW analysis, no time-lag testing
    if nargin<2
        suffix = 'pw';
    end
end
method = containers.Map(1,'RA');
method(2) = 'CO';
method(3) = 'GC';
method(4) = 'MI';
method(5) = 'TE';
method(6) = 'CCM';
dream = strcat('data/data_dream_', num2str(num_nodes),'_');
data_pre = strcat('data/result_dream_', num2str(num_nodes),'_');
if strcmp(lagtest,'')
    data_suf = strcat('_correct_', suffix);
else
    data_suf = strcat('_correct_', suffix, 'lagtest', lagtest);
end
out = [];
if strcmp(suffix, 'pw')
    for i=1:5
        dream_nam = strcat(dream, num2str(i), '_correct');
        data_nam = strcat(data_pre, num2str(i), data_suf);
        load(dream_nam);
        load(data_nam);
        original = zeros(num_nodes);
        for j=1:size(grn,1)
            original(grn(j,1),grn(j,2)) = 1;
        end
        found = pw_out{6};
%         in_outdegs = [];
%         in_indegs = [];
%         out_outdegs = [];
%         out_indegs = [];
        for k=1:length(found) %"random" also included
            detected = found{k};
            in_outdeg = [];
            in_indeg = [];
            out_outdeg = [];
            out_indeg = [];
            for j=1:size(grn,1)
                x = detected(j,1);
                y = detected(j,2);
                in_outdeg = cat(1,in_outdeg,sum(original(x,:)));
                in_indeg = cat(1,in_indeg,sum(original(:,x)));
                out_outdeg = cat(1,out_outdeg,sum(original(y,:)));
                out_indeg = cat(1,out_indeg,sum(original(:,y)));
            end
%             in_outdegs = cat(2,in_outdegs,in_outdeg);
%             in_indegs = cat(2,in_indegs,in_indeg);
%             out_outdegs = cat(2,out_outdegs,out_outdeg);
%             out_indegs = cat(2,out_indegs,out_indeg);
            [counts, bins] = hist(cat(2,in_outdeg,in_indeg,out_outdeg,out_indeg));
            plot(bins, counts);
            title(strcat('Frequency Distributions of Node-Degrees of Detected Links for', method(k), ' on DREAM4-', num2str(num_nodes), ' #', num2str(i)));
            xlabel('Degree');
            ylabel('Distribution');
            legend('in-outdeg','in-indeg','out-outdeg','out-indeg');
            print(['../plots/plots_analysis/dream_' num2str(num_nodes) '_' num2str(i),'_correct_' method(k) 'linkfreq'], '-dpng');
        end
    end
end

