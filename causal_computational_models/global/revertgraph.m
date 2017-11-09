function [ nodematrix ] = revertgraph( linkmatrix )
%REVERTGRAPH Reverts a link-graph into its node-graph
n = sqrt(size(linkmatrix,1));
ranks = pagerank(linkmatrix);%stationary_distribution_ctmc(linkmatrix);%
nodematrix = zeros(n);
for i=1:n
    for j=1:n
        if (i==j)
            continue
        end
        nodematrix(i,j) = ranks(sub2ind([n,n],i,j));
    end
end
%nodematrix
end

