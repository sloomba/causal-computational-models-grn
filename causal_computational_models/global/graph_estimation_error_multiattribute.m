function [ error_tot ] = graph_estimation_error_multiattribute( rho, data_new, theta_and, lambda )
%GRAPH_ESTIMATION_ERROR Error function for graph_estimation.m
if nargin<4
    nargin
    lambda = 0; %regulariser
end
dig_lap = diag(sum(theta_and,2)) - theta_and;
a = size(data_new,3);
n = size(data_new,1);
error_tot = 0;
for i=1:a
    error = (data_new(:,:,i) - rho(i,1)*expm(-rho(i,2)*dig_lap)).^2;
    error(1:n+1:n*n) = 0;
    error_tot = error_tot+sum(sum(error));
end
error_tot = error_tot+lambda*sum(sum(theta_and));
end