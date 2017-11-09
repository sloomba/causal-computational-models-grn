function [ clusters, t_f ] = spectral_clustering_links( expn, grn, alpha, technique, num, range )
%SPECTRAL_CLUSTERING
if nargin<6
    range = 10; %less than time-series length
    if nargin<5
        num = 4; %less than 8
        if nargin<4
            technique = 'Laplacian';
            if nargin<3
                alpha = 0.01;
            end
        end
    end
end
s = []; %link matrix
for i=1:size(expn,1)
    for j=1:i-1
        co1 = xcorr(expn(i,:),expn(j,:));
        co1 = co1(1:floor((length(co1)-1)/2));
        co2 = xcorr(expn(j,:),expn(i,:));
        co2 = co2(1:floor((length(co2)-1)/2));
        s = cat(1,s,co1);
        s = cat(1,s,co2);
    end
end
s = s(:,1:range);
%s = compute_mapping(s,technique,range);
sim = pdist(s); %distance matrix
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
clusters = 2*kmeans(eigvec,num);
s2 = compute_mapping(s,technique);
index = 1;
for i=1:size(expn,1)
    for j=1:i-1
        flag1 = 0;
        flag2 = 0;
        for k=1:size(grn,1)
            if all(grn(k,1:2)==[i,j])
                flag1 = 1;
            end
            if all(grn(k,1:2)==[j,i])
                flag2 = 1;
            end
            if flag1 && flag2
                break;
            end
        end
        clusters(index) = clusters(index)-flag1;
        clusters(index+1) = clusters(index+1)-flag2;
        index = index+2;
    end
end
gscatter(s2(:,1),s2(:,2),clusters,'bbrrggmmcckkyy','o*o*o*o*o*o*o*');
t_f = [];
for i=1:num
    t = length(find(clusters==2*i-1));
    f = length(find(clusters==2*i));
    t_f = cat(1,t_f,[100*t/(t+f),100*f/(t+f)]);
end

