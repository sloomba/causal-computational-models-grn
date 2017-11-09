function [ mintheta, theta_ands, errors, minrho, minidx, minval ] = graph_estimation_multiattribute( data, lambda, limit, max_iters )
%GRAPH_ESTIMATION Intrinsic graph estimation method using digraph
%Laplacian. Data is an nXnXa data matrix, where 'n' is number of nodes,
%'a' is number of attributes. 'mintheta' is the estimated nXn adjacency
%matrix (directed graph).
%SEE NODA ET AL. We extend it to multiple attributes.
warning('off','all');
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
a = size(data,3);
while true
    data_new = data + repmat(r*diag(ones(1,n)),[1,1,a]);
    flag = true;
    for i=1:a
        if det(data_new(:,:,i))== 0
            r = r+1;%abs(normrnd(1,0.25));
            flag = false;
            break;
        end
    end
    if flag
        break;
    end
end
display(strcat('Data matrix initialised with r=',num2str(r)))
minidx = 1;
minval = Inf;
minrho = [];
mintheta = [];
theta_ands = [];
errors = [];
lambda
limit
for m=1:n*(n-1)
    m
    for iter=1:max_iters
        temps = [];
        thetas = [];
        for i=1:a
            theta = abs(logm(data_new(:,:,i)));
            theta(1:n+1:n*n) = 0;
            temps = cat(1,temps,sort(reshape(theta,1,n*n),'descend'));
            thetas = cat(3,thetas,theta);
        end
        tryidx = m;
        theta_and = ones(n);
        while true
            for i=1:a
                val = temps(i,tryidx);
                theta = thetas(:,:,i);
                theta(theta<val) = 0;
                theta(theta>=val) = 1; %find current theta
                theta_and = theta_and & theta;
            end
            if sum(sum(theta_and))>=m
                break;
            else
                tryidx = tryidx+1;
                theta_and = ones(n);
            end
        end
        %display(strcat('Margin of disagreement between attributes for m=',num2str(m),' is del=',num2str(tryidx-m))); %uncommment if verbose
        rho = fminsearch(@(rho) graph_estimation_error_multiattribute(rho,data_new,theta_and,lambda),rand(a,2)); %argmin over rho
        dig_lap = diag(sum(theta_and,2)) - theta_and;
        update = [];
        for i=1:a
            update = cat(3,update,diag(diag(rho(i,1)*expm(-rho(i,2)*dig_lap))));
        end
        data_new = data + update; %update data matrix
    end
    error = graph_estimation_error_multiattribute(rho,data_new,theta_and,lambda);
    if (error<minval)
        minidx = m;
        minval = error;
        minrho = rho;
        mintheta = theta_and;
    end
    theta_ands = cat(3,theta_ands,theta_and);
    errors = cat(2,errors,error);
end
warning('on','all');
end

