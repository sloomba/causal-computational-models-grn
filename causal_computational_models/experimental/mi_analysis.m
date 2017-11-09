function [ out ] = mi_analysis( expn_data, grn, num_nodes )
%Mutual Information analysis for given expression data and ground truth GRN.
all_mis = [];
mis_links = [];
all_links = containers.Map(-1,containers.Map(-1,0));
all_links.remove(-1);
for i=1:num_nodes
    temp = containers.Map(1,0);
    for j=2:num_nodes
        if i==j
            continue;
        else
            temp(j) = 0;
        end
    end
    all_links(i) = temp;
end
for link=1:length(grn)
    mis = [];
    y = grn(link,1)
    x = grn(link,2)
    temp = all_links(y);
    temp.remove(x);
    all_links(y) = temp;
    for i=0:15
        i
        mis = cat(1,mis,mutual_information(expn_data(:,y),expn_data(:,x),i));
    end
    all_mis = cat(2,all_mis,mis);
    mis_links = cat(2,mis_links,[y;x]);
end
all_non_mis = [];
non_mis_links = [];
keylist = keys(all_links);
for i=1:length(keylist)
    y = keylist{i}
    temp_dict = all_links(y);
    temp = keys(temp_dict);
    for j=1:length(temp)
        x = temp{j}
        non_mis = [];
        for k=0:15
            k
            non_mis = cat(1,non_mis,mutual_information(expn_data(:,y),expn_data(:,x),k));
        end
        all_non_mis = cat(2,all_non_mis,non_mis);
        non_mis_links = cat(2,non_mis_links,[y;x]);
    end
end
out = {all_mis, mis_links, all_non_mis, non_mis_links};
end
