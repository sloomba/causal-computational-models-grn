function [ out ] = te_analysis( expn_data, grn, num_nodes )
%Transfer Entropy analysis for given expression data and ground truth GRN.
all_tes = [];
tes_links = [];
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
    tes = [];
    y = grn(link,1)
    x = grn(link,2)
    temp = all_links(y);
    temp.remove(x);
    all_links(y) = temp;
    for i=1:16
        i
        tes = cat(1,tes,transfer_entropy(expn_data(:,y),expn_data(:,x),i));
    end
    all_tes = cat(2,all_tes,tes);
    tes_links = cat(2,tes_links,[y;x]);
end
all_non_tes = [];
non_tes_links = [];
keylist = keys(all_links);
for i=1:length(keylist)
    y = keylist{i}
    temp_dict = all_links(y);
    temp = keys(temp_dict);
    for j=1:length(temp)
        x = temp{j}
        non_tes = [];
        for k=1:16
            k
            non_tes = cat(1,non_tes,transfer_entropy(expn_data(:,y),expn_data(:,x),k));
        end
        all_non_tes = cat(2,all_non_tes,non_tes);
        non_tes_links = cat(2,non_tes_links,[y;x]);
    end
end
out = {all_tes, tes_links, all_non_tes, non_tes_links};
end

