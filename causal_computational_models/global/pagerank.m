function [ ranks ] = pagerank( matrix, d, epsilon, alpha )
%PAGERANK matrix could be an adjacency matrix or a weighted matrix
if nargin<4
    alpha = 0.00000001;
    if nargin<3
        epsilon = 0.0001;
        if nargin<2
            d = 0.85;
        end
    end
end
n = size(matrix,1);
if all(all(matrix==0))
    matrix(matrix==0) = alpha; %if a zero matrix, correct it
else
    maxval = max(max(matrix));
    newmat = matrix;
    newmat(newmat==0) = maxval;
    minval = min(min(newmat));
    matrix(matrix==0) = alpha*minval;
end
weights = repmat(sum(matrix,2),1,n);
ranks = 2*ones(1,n); %initialise ranks to 1
old_ranks = ones(1,n);
while all(abs(ranks-old_ranks)>epsilon)
    old_ranks = ranks;
    ranks = 1-d + d*old_ranks*(matrix./weights);
end
end

