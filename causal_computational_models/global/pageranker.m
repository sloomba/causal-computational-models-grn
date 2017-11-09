function [ nodematrix, all_adjs, all_errors ] = pageranker( pairwisematrix, epsilon, maxiters, eta )
%PAGERANKER Iterative pagerank for graph construction
if nargin<4
    eta = 0.5;
    if nargin<3
        maxiters = 2;
        if nargin<2
            epsilon = 0.0001;
        end
    end
end
n = size(pairwisematrix,1);
all_adjs = [];
all_errors = [];
iters = 0;
curr_pairwisematrix = pairwisematrix;
temp = pairwisematrix;
maxv = max(max(temp));
temp(temp==0)=maxv;
minv = min(min(temp));
if any(any(pairwisematrix<0))
    display('HAWWWWW');
    return
end
error =  0;
while iters<maxiters
    last_pairwisematrix = curr_pairwisematrix;
    curr_pairwisematrix = (curr_pairwisematrix.*revertgraph(invertgraph(curr_pairwisematrix)));
    %curr_pairwisematrix = revertgraph(invertgraph(curr_pairwisematrix));
    %curr_pairwisematrix = (1-eta)*curr_pairwisematrix + eta*revertgraph(invertgraph(curr_pairwisematrix));
%     temp = curr_pairwisematrix;
%     new_maxv = max(max(temp));
%     temp(temp==0)=new_maxv;
%     new_minv = min(min(temp));
%     curr_pairwisematrix = ((maxv-minv)/(new_maxv-new_minv))*curr_pairwisematrix + ((minv*new_maxv-maxv*new_minv)/(new_maxv-new_minv))*ones(n);
    error = norm(curr_pairwisematrix-last_pairwisematrix,'fro')/n;
    %inner_error = cat(1,inner_error,error);
    if error<epsilon
        break
    end
    iters = iters+1;
end
iters
topm = sort(reshape(curr_pairwisematrix,1,n*n),'descend');
for m=1:n*(n-1)
    topmv = topm(m);
    nodematrix = curr_pairwisematrix;
    nodematrix(nodematrix<topmv) = 0;
    nodematrix(nodematrix~=0) = 1;
    all_adjs = cat(3,all_adjs,nodematrix);
    all_errors = cat(1,all_errors,error);
end
idx = all_errors==min(all_errors);
nodematrix = all_adjs(:,:,idx);
nodematrix = nodematrix(:,:,1);
end

