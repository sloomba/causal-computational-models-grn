function [ edges, non_edges ] = dim_analysis( expn, grn, technique, dims )
%does dimensionality reduction analysis
if nargin<4
    dims = 2;
end
new_expn = [];
for i=1:size(expn,1)
    for j=1:i-1
        new_expn = cat(1,new_expn,xcorr(expn(i,:),expn(j,:)));
    end
end
new_expn_dim = compute_mapping(new_expn,technique,dims);
index = 1;
edges = [];
non_edges = [];
for i=1:size(expn,1)
    for j=1:i-1
        flag = 0;
        for k=1:size(grn,1)
            if all(grn(k,1:2)==[i,j]) || all(grn(k,1:2)==[j,i])
                flag = 1;
                break;
            end
        end
        if flag
            edges = cat(1,edges,new_expn_dim(index,:));
        else
            non_edges = cat(1,non_edges,new_expn_dim(index,:));
        end
        index = index+1;
    end
end
plot(edges(:,1),edges(:,2),'ob');
hold all;
plot(non_edges(:,1),non_edges(:,2),'*r');
end

