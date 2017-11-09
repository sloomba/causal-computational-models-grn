function [ clusters ] = spectral_clustering_nodes( expn, grn, alpha, technique, num )
%SPECTRAL_CLUSTERING
if nargin<5
    num = 7; %less than 8
    if nargin<4
        technique = 'Laplacian';
        if nargin<3
            alpha = 0.01;
        end
    end
end
sim = pdist(expn); %distance matrix
sim = squareform(sim);
aff = exp(-sim*alpha); %affinity matrix
n = size(aff,1);
aff(1:n+1:n^2) = 0;
diavec = sum(aff,2);
dia = zeros(n); %diagonal matrix
dia(1:n+1:n^2) = diavec(1:n);
invdia = sqrt(inv(dia));
lap = invdia*aff*invdia; %laplacian matrix
[eigvec,~] = eig(lap);
eigvec = eigvec(:,1:num);
for i=1:n
    eigvec(i,:) = eigvec(i,:)/sqrt(sum(eigvec(i,:).*eigvec(i,:)));
end
%clusters = kmeans(eigvec,num);
clusters = kmeans(expn,num);
s = compute_mapping(expn,technique);
subplot(1,2,1);
gscatter(s(:,1),s(:,2),clusters,'brgmcky','o');
xlabel('Dim 1');
ylabel('Dim 2');
for i=1:size(expn,1)
    text(s(i,1)+0.005,s(i,2)+0.005,num2str(i));
end
subplot(1,2,2);
scatter(grn(:,1),grn(:,2),'b','o');
xlabel('From node');
ylabel('To node');
end
    
