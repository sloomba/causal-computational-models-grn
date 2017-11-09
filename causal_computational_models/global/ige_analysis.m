function [ out ] = ige_analysis( grn, exp, metric, figname )
%IGE_ANALYSIS Intrinsic Graph Estimation Analysis
addpath('../helpers');
n = size(exp,1);
len = size(exp,2);
if nargin<3
    metric = zeros(n); %if metric matrix is not give, only then is a "corrleogram" matrix created and used as metric
    for i=1:n
        for j=1:i-1
            dummy = abs(xcov(exp(i,:),exp(j,:)));
            a = max(dummy(1:len));
            b = max(dummy(len:end));
            metric(j,i) = a;
            metric(i,j) = b;
        end
    end
end
adj = edgelist_to_adjmat(grn,n);
size(metric)
[adj_found, all_adjs, all_errors] = graph_estimation(metric); %juggler(metric); %can also append '_multiattribute'
adj_found
adj_found_list = adjmat_to_edgelist(adj_found);
roc = [];
precision = [];
recall = [];
fscore = [];
for i=1:size(all_adjs, 3)
    curr_adj = all_adjs(:,:,i);
    tp = sum(sum(adj&curr_adj));
    fp = sum(sum(curr_adj-adj&curr_adj));
    tn = sum(sum(~(adj|curr_adj)))-n;
    fn = sum(sum(~curr_adj-~(adj|curr_adj)));
    roc = cat(1,roc,[fp/(fp+tn),tp/(tp+fn)]);
    precision = cat(1,precision,tp/(tp+fp));
    recall = cat(1,recall,tp/(tp+fn));
    fscore = cat(1,fscore,2*precision(end)*recall(end)/(precision(end)+recall(end)));
end
[recall, sort_order] = sort(recall);
precision = precision(sort_order);
fscore = fscore(sort_order);
[~,idxs]=sort(roc(:,2));
roc = roc(idxs,:);
i = 1;
last = roc(1,1);
while i<size(roc,1)
    if roc(i,1)<last
        roc(i,:) = [];
    else
        last = roc(i,1);
        i = i+1;
    end
end
auroc = trapz(roc(:,1),roc(:,2));
tp = sum(sum(adj&adj_found));
fp = sum(sum(adj_found-adj&adj_found));
tn = sum(sum(~(adj|adj_found)))-n;
fn = sum(sum(~adj_found-~(adj|adj_found)));
p = tp/(tp+fp);
r = tp/(tp+fn);
f = 2*p*r/(p+r);
a = (tp+tn)/(tp+tn+fp+fn);
out = {p, r, f, a, adj_found, roc, auroc, precision, recall, fscore, all_adjs, all_errors};
if nargin>=4
    figure(2); clf;
    subplot(1,2,1);
    hold;
    plot(precision, 'color', 'r');
    plot(recall, 'color', 'g');
    plot(fscore, 'color', 'b');
    plot(p*ones(size(precision)), 'color', 'r');
    text(length(precision), p+0.02, num2str(p));
    plot(r*ones(size(recall)), 'color', 'g');
    text(length(recall), r+0.02, num2str(r));
    plot(f*ones(size(fscore)), 'color', 'b');
    text(length(fscore), f+0.02, num2str(f));
    legend('Precision', 'Recall', 'f-score');
    xlabel('Links');
    ylabel('Score');
    hold;
    title(strcat('Precision, Recall and f-score plot for IGE with given metric'));
    subplot(1,2,2);
    plot(roc(:,1),roc(:,2));
    xlabel('False Positive Rate');
    ylabel('True Positive Rate');
    text(0.5, 0.2, num2str(auroc));
    title('ROC plot');
    print(strcat(figname,'ROCplotIGE'),'-dpng');
    figure(2); clf;
    hold;
    scatter(grn(:,1),grn(:,2),'bo');
    scatter(adj_found_list(:,1),adj_found_list(:,2),'g*');
    legend('Ground', 'Found');
    xlabel('From node');
    ylabel('To node');
    hold;
    title('Scatter plots for given metric');
    print(strcat(figname,'ScatterIGE'),'-dpng');
end
end

