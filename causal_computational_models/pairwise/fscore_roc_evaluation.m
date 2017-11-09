function [ precision, recall, fscore, fpr, tpr ] = fscore_roc_evaluation( num_nodes, ground, found )
%Evaluates the discovered causalities.
tp = 0;
fp = 0;
num_links = size(found,1);
precision = zeros(num_links,1);
recall = zeros(num_links,1);
fscore = zeros(num_links,1);
fpr = zeros(num_links,1);
tpr = zeros(num_links,1);
for i=1:size(found,1)
    flag = true;
    for j=1:size(ground,1)
        if all(ground(j,1:2)==found(i,1:2))
            tp = tp+1;
            flag = false;
            break;
        end
    end
    if flag
        fp = fp+1;
    end
    fn = size(ground,1)-tp;
    tn = num_nodes*(num_nodes-1)-size(ground,1)-fp;
    precision(i) = tp/(tp+fp);
    recall(i) = tp/(tp+fn);
    fscore(i) = 2*tp/(2*tp+fp+fn);
    fpr(i) = fp/(fp+tn);
    tpr(i) = tp/(tp+fn);
end
end

