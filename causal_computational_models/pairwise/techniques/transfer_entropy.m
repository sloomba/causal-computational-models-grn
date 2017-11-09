function [ trans_ent ] = transfer_entropy( y, x, lag )
%Computes transfer entropy y --> x with given lag parameters.
assert(length(x)==length(y));
symbolsy = unique(y);
symbolsx = unique(x);
sy = containers.Map(symbolsy,1:length(symbolsy));
sx = containers.Map(symbolsx,1:length(symbolsx));
Px = ones(length(symbolsx),1); %if no Laplace Smoothing, zeros
Pxy = ones(length(symbolsx),length(symbolsy));
Pxx = ones(length(symbolsx),length(symbolsx));
Pxxy = ones(length(symbolsx),length(symbolsx),length(symbolsy));
count = idivide(length(x), int32(21));
ym = y;
xm = x;
for k=1:count
    x = xm((k-1)*21 + 1: k*21);
    y = ym((k-1)*21 + 1: k*21);
    for i=1:length(x)
        Px(sx(x(i))) = Px(sx(x(i)))+1;
        Pxy(sx(x(i)),sy(y(i))) = Pxy(sx(x(i)),sy(y(i)))+1;
        if i<=length(x)-lag
            Pxx(sx(x(i+lag)),sx(x(i))) = Pxx(sx(x(i+lag)),sx(x(i)))+1;
            Pxxy(sx(x(i+lag)),sx(x(i+lag-1)),sy(y(i))) = Pxxy(sx(x(i+lag)),sx(x(i+lag-1)),sy(y(i)))+1;
        end
    end
end
countx = sum(Px);
countxy = sum(sum(Pxy));
countxx = sum(sum(Pxx));
countxxy = sum(sum(sum(Pxxy)));
trans_ent = 0;
for i=1:length(symbolsx)
    for j=1:length(symbolsx)
        vx = Px(j)/countx;
        vxx = Pxx(i,j)/countxx;
        for k=1:length(symbolsy)
            vxxy = Pxxy(i,j,k)/countxxy;
            vxy = Pxy(j,k)/countxy;
            to_add = vxxy*log2((vxxy*vx)/(vxy*vxx));
            if ~isnan(to_add)
                trans_ent = trans_ent + to_add;
            end
        end
    end
end
if trans_ent<0 %Transfer Entropy can be negative!
    trans_ent = 0;
end
end