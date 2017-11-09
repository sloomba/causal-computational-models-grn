function [ f_score, precision, recall, fpr ] = gc_analysis( data, ground_truth, alpha, max_lag, str_idx, stp_idx )
%Given the gene expression level data (column-wise), this function
%evaluates causality between signals using Granger Causality, and checks
%the obtained network against the ground truth, thus return the f-score
%values.
tic
causalities = [];
correlation = [];
for i=max(1,str_idx):min(size(data,2),stp_idx)
    %i
    for j=max(1,str_idx):min(size(data,2),stp_idx)
        if (i==j)
            continue;
        end
        [F, c_v] = granger_cause(data(:,i), data(:,j), alpha, max_lag);
        if (F>c_v)
            causalities = cat(1, causalities, [j, i]);
	    correlation = cat(1, correlation, corr2(data(:, i), data(:, j)));
        end
    end
end
found = zeros(1,size(ground_truth,1));
correct_gt = zeros(1,size(ground_truth,1));
correct_nu = zeros(1,size(causalities,1));
for i=1:size(ground_truth,1)
    if (all(ground_truth(i,1:2)>=str_idx) && all(ground_truth(i,1:2)<=stp_idx))
        found(i) = 1;
        for j=1:size(causalities,1)
            if (correct_nu(j)==0)
                if (all(ground_truth(i,1:2)==causalities(j,:)))
                    correct_nu(j)=1;
                    correct_gt(i)=1;
                    break;
                end
            end
        end
    end
end
%dlmwrite(strcat('causalities_', num2str(size(data, 2)), '_nodes_lag_', num2str(max_lag), '_aplha_', num2str(alpha), '.tsv'), causalities, 'delimiter', '\t');
tp = sum(correct_gt);
precision = tp/size(causalities,1);
recall = tp/sum(found);
f_score = 2*precision*recall/(precision+recall);
n = min(size(data,2),stp_idx) - max(1,str_idx) + 1;
tn = n*(n-1) - sum(found) + tp - size(causalities,1);
fpr = (size(causalities,1)-tp)/(tn+(size(causalities,1)-tp));
plot(correlation);
print(strcat('results/', num2str(n), '_nodes/correlation_', num2str(max_lag), '_maxlag/', 'alpha_', num2str(alpha), '.png'), '-dpng');
toc

% n*n - n - found + tp - causalities = tn
% tp + fp = causalities
% found = tp + fn
% required = tn

