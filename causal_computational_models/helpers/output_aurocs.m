function [ out ] = output_aurocs( num_nodes, suffix, lagtest )
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
data_pre = strcat('../data/result_dream_', num2str(num_nodes),'_');
if strcmp(lagtest,'')
    data_suf = strcat('_correct_', suffix);
else
    data_suf = strcat('_correct_', suffix, 'lagtest', lagtest);
end
out = [];
if strcmp(suffix, 'pw')
    for i=1:5
        data_nam = strcat(data_pre, num2str(i), data_suf);
        load(data_nam);
        fpr = pw_out{4};
        tpr = pw_out{5};
        curr = [];
        for j=1:length(fpr) %"random" also included
            curr = cat(1, curr, abs(trapz(fpr{j}, tpr{j})));
        end
        out = cat(2, out, curr);
    end
else
    for i=1:5
        data_nam = strcat(data_pre, num2str(i), data_suf);
        load(data_nam);
        curr = [];
        for j=2:6 %"random" not included
            currval = ige_out{j};
            curr = cat(1, curr, currval{7});
        end
        out = cat(2, out, curr);
    end
end

