function [ nodematrix, all_adjs, all_errors ] = juggler( pairwisematrix, epsilon, maxiters, eta )
%JUGGLER Iterative pagerank for graph construction
if nargin<4
    eta = 0.5;
    if nargin<3
        maxiters = 20;
        if nargin<2
            epsilon = 0.0001;
        end
    end
end
n = size(pairwisematrix,1);
all_adjs = zeros(n, n, n*(n-1));
all_errors = zeros(n*(n-1), 1);
% all_adjs = [];
% all_errors = [];
%inner_errors = [];
parfor m=1:n*(n-1)
    m
    iters = 1;
    curr_pairwisematrix = pairwisematrix;
    %inner_error = [];
    while iters<maxiters
        nodematrix = revertgraph(invertgraph(curr_pairwisematrix));
        topm = sort(reshape(nodematrix,1,n*n),'descend');
        topm = topm(m);
        nodematrix = (1 + tanh(nodematrix-topm*ones(n)))/2;
        last_pairwisematrix = curr_pairwisematrix;
        curr_pairwisematrix = (1-eta)*curr_pairwisematrix + eta*revertgraph(invertgraph(nodematrix));
        error = norm(curr_pairwisematrix-last_pairwisematrix,'fro')/n
        %inner_error = cat(1,inner_error,error);
        if error<epsilon
            break
        end
        iters = iters+1;
    end
    %inner_errors = cat(2,inner_errors,inner_error);
    iters
    topm = sort(reshape(nodematrix,1,n*n),'descend');
    topm = topm(m);
    nodematrix(nodematrix<topm) = 0;
    nodematrix(nodematrix~=0) = 1;
    % all_adjs = cat(3,all_adjs,nodematrix);
    all_adjs(:, :, m) = nodematrix;
    % all_errors = cat(1,all_errors,error);
    all_errors(m) = error;
end
idx = all_errors==min(all_errors);
nodematrix = all_adjs(:,:,idx);
end

