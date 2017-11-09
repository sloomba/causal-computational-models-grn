function [ out ] = edgelist_to_adjmat( in, n )
%EDGELIST_TO_ADJMAT 'n' is number of nodes.
out = zeros(n);
for i=1:size(in,1)
    out(in(i,1),in(i,2)) = 1;
end
end

