function [ mut_info ] = mutual_information( y, x , lag )
%Computes mutual information between y and (-->) x, with given lag parameters.
assert(length(x)==length(y));
symbolsy = unique(y);
symbolsx = unique(x);
sy = containers.Map(symbolsy,1:length(symbolsy));
sx = containers.Map(symbolsx,1:length(symbolsx));
Px = ones(length(symbolsx),1); %zeros, if no Laplace Smoothing
Py = ones(length(symbolsy),1);
Pxy = ones(length(symbolsx),length(symbolsy));
xm = x;
ym = y;
count = double(idivide(length(x), int32(21)));
for k=1:count
    x = xm((k-1)*21 + 1: k*21);
    y = ym((k-1)*21 + 1: k*21);
    if lag>0    %y-->x
        x = x(lag+1:end);
        y = y(1:end-lag);
    else        %x-->y
        y = y(lag+1:end);
        x = x(1:end-lag);
    end
    for i=1:length(x)
        Px(sx(x(i))) = Px(sx(x(i)))+1;
        Py(sy(y(i))) = Py(sy(y(i)))+1;
        Pxy(sx(x(i)),sy(y(i))) = Pxy(sx(x(i)),sy(y(i)))+1;
    end
end
countx = sum(Px);
county = sum(Py);
countxy = sum(sum(Pxy));
mut_info = 0;
for i=1:length(symbolsx)
    for j=1:length(symbolsy)
        to_add = Pxy(i,j)*log2((Pxy(i,j)*countx*county)/(Px(i)*Py(j)*countxy))/countxy;
        if ~isnan(to_add)
            mut_info = mut_info + to_add;
        end
    end
end
end