function [ mintheta, thetas, errors, minrho, minidx, minval ] = graph_estimation_new( data, lambda, limit, max_iters )
%GRAPH_ESTIMATION Intrinsic graph estimation method using digraph
%Laplacian. Data is an nXn data matrix, where 'n' is number of nodes.
%'mintheta' is the estimated nXn adjacency matrix (directed graph).
%warning('off','all');
%SEE NODA ET AL. We introduce soft thresholding which improves results.
if nargin<4
    max_iters = 4;
    if nargin<3
        limit = Inf;
        if nargin<2
            lambda = 0;
        end
    end
end
r = 1;
n = size(data,1);
while true
    data_new_ori = data + r*diag(ones(1,n));
    if det(data_new_ori)~=0
        break;
    else
        r = r+1;%abs(normrnd(1,0.25));
    end
end
display(strcat('Data matrix initialised with r=',num2str(r)));
rhos = [];
thetas = [];
errors = [];
lambda
limit
options = optimset('MaxFunEvals',100000,'MaxIter',100000);
parfor m=1:min(n*(n-1),limit)
    m
    data_new = data_new_ori;
    for iter=1:max_iters
        theta = abs(logm(data_new));
        theta(1:n+1:n*n) = 0;
        temp = sort(reshape(theta,1,n*n),'descend');
        val = temp(m);
        theta = (1 + tanh(theta-val*ones(n)))/2; %comment for hard thresholding
        %theta(theta<val) = 0; %uncomment block for hard thresholding
        %theta(theta>=val) = 1; %uncomment block for hard thresholding %find current theta
        rho = fminsearch(@(rho) graph_estimation_error(rho,data_new,theta,lambda),[rand,rand],options); %argmin over rho
        dig_lap = diag(sum(theta,2)) - theta;
        data_new = data + diag(diag(rho(1)*expm(-rho(2)*dig_lap))); %update data matrix
    end
    error = graph_estimation_error(rho,data_new,theta,lambda);
    rhos(m,:) = rho;
    theta(theta<0.5) = 0; %comment block for hard thresholding
    theta(theta>=0.5) = 1; %comment block for hard thresholding
    thetas(:, :, m) = theta;   % cat(3,thetas,theta);
    errors(m) = error;% cat(1,errors,error);
end
[minval, minidx] = min(errors);
mintheta = thetas(:,:,minidx);
minrho = rhos(minidx,:);
%warning('on','all');
end