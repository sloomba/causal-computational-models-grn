function [ out ] = metriclist_to_metricmat( list )
%METRICLIST_TO_METRICMAT
n = size(list,1);
i = 1;
while i<n
    if n==i*(i-1)
        break;
    end
    i = i+1;
end
out = zeros(i);
for j=1:n
    out(list(j,2),list(j,3)) = list(j,1);
end
end

