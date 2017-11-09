function [ linkmatrix ] = invertgraph( nodematrix )
%INVERTGRAPH Inverts a node-graph into its link-graph
n = size(nodematrix,1);
ranks = pagerank(nodematrix);%stationary_distribution_ctmc(nodematrix);%
linkmatrix = zeros(n*n);
for i=1:n
    for j=1:n
        if (i==j)
            continue
        end
        for k=1:n
            if (i==k)
                continue
            end
            from = sub2ind([n,n],j,i);
            to = sub2ind([n,n],i,k);
            linkmatrix(from,to) = ranks(i)*(nodematrix(j,i)*nodematrix(i,k));
        end
    end
end
%linkmatrix
end

