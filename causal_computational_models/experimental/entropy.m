function [ ent ] = entropy( signal, bin_size )
%Calculates entropy of a given (exact or binned) signal.
ent = 0;
counts = containers.Map(-1,-1);
counts.remove(-1);
total = 0;
for i=1:length(signal)
    bin = floor((signal(i))/bin_size);
    if counts.isKey(bin)
        counts(bin) = counts(bin)+1;
        total = total+1;
    else
        counts(bin) = 1;
    end
end
keylist = keys(counts);
for i=1:length(keylist)
    prob = counts(keylist{i})/total;
    if prob>0
        ent = ent - prob*log2(prob);
    end
end
end

