function [ p ] = ptm_convergence( p, k, epsilon )
%PTM_CONVERGENCE
%THIS SCRIPT WAS DEVELOPED IN THE EARLIER STAGES OF THE PROJECT. CAN BE IGNORED.
if nargin<3
    epsilon = 0.01;
end
n=size(p,1);
oldp = zeros(n);
while true
    p = p/(det(p)^(1/n)) %ensure determinant=1
    for i=1:n
        p(i,i)=0;
    end
    error = sum(sum(abs(p-oldp)))/n
    if error <epsilon
        break
    end
    oldp = p;
    p = p*p;
end
q = reshape(p,1,n*n);
q = sort(q,'descend');
topk = q(k);
p(p<topk) = 0;
p(p~=0) = 1;
end

