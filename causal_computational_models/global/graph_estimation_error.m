function [ error ] = graph_estimation_error( rho, data_new, theta, lambda )
%GRAPH_ESTIMATION_ERROR Error function for graph_estimation.m
if nargin<4
    nargin
    lambda = 0; %regulariser
end
dig_lap = diag(sum(theta,2)) - theta;
error = (data_new - rho(1)*expm(-rho(2)*dig_lap)).^2;
n = size(theta,1);
error(1:n+1:n*n) = 0;
error = sum(sum(error))+lambda*sum(sum(theta));
end