function [ x ] = stationary_distribution_ctmc( q )
%STATIONARY_DISTRIBUTION_CTMC Finds the stationary distribution of a
%continuous time Markov chain, given the rate matrix q.
%THIS SCRIPT WAS DEVELOPED IN THE EARLIER STAGES OF THE PROJECT. CAN BE IGNORED.
n = size(q,2);
for i=1:n
    q(i,i) = q(i,i)-sum(q(i,:));
end
a = cat(1,q',ones(1,n));
b = cat(1,zeros(n,1),1);
x = transpose(linsolve(a,b));
end

