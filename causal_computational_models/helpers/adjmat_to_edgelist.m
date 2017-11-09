function [ out ] = adjmat_to_edgelist( in )
%ADJMAT_TO_EDGELIST
n = size(in,1);
out = [];
for i=1:n
    for j=1:n
        if i==j
            continue;
        else
            if in(i,j)~=0
                out = cat(1,out,[i,j,in(i,j)]);
            end
        end
    end
end
end

