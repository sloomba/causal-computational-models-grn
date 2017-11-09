function [ precision, recall, fscore, fpr, tpr, sorted_all_values, all_values ] = roctest( exp, grn, technique, max_lag, figname)
%Does ROC analysis using given pairwise technique, assuming
%num_nodes > num_time_steps.
%Finds optimal lag using AIC/BIC
%Doesn't count for self regulation
%techniques = 'ra', 'co', 'gc', 'mi', 'te', 'cm'
if nargin<5
    figname = '';
    if nargin<4
        max_lag = 10;
        if nargin<3
            technique = 'ra';
        end
    end
end
technique
warning off;
% if (size(exp,1)>=size(exp,2))
%     exp = exp';
% end
if strcmp(technique, 'mi') || strcmp(technique, 'te')
    exp = cont_to_disc(exp, 20);
end
num_nodes = size(exp,1);
all_values = [];
if strcmp(technique,'ra')
    for i=1:num_nodes
        for j=1:num_nodes
            if i==j
                continue;
            else
                all_values = cat(1,all_values,[i,j]);
            end
        end
    end
    all_values = all_values(randperm(size(all_values,1)),:);
    sorted_all_values = all_values;
elseif strcmp(technique,'cm')
    for i=1:num_nodes
        for j=1:i-1
            i
            j
            c1s = 0;
            c2s = 0;
            [~, ~ , ~ , SugiY , SugiX , origY , origX] = convergent_cross_map(exp(i,:)',exp(j,:)',1,2);
            SugiX(~isfinite(SugiX)) = 0;
            SugiY(~isfinite(SugiY)) = 0;
            c1 = corrcoef(SugiX, origX);
            if (c1(1, 2)>0)
                c1s = c1(1, 2);
            end
            c2 = corrcoef(SugiY, origY);
            if (c2(1, 2)>0)
                c2s = c2(1, 2);
            end
            all_values = cat(1,all_values,[c1s,i,j,1]);
            all_values = cat(1,all_values,[c2s,j,i,1]);
            if ~strcmp(figname,'')
                plot(origX, SugiX, 'or')
                hold all
                plot(origY, SugiY, 'og')
                line = [min(min(min(origX),min(SugiX)),min(min(origY),min(SugiY))), max(max(max(origX),max(SugiX)),max(max(origY),max(SugiY)))];
                plot(line,line,'b')
                hold off
                legend('X', 'Y');
                xlabel('Original Signal');
                ylabel('Reconstructed Signal');
                title(strcat('Scatter plots for reconstructed signals using CCM (', num2str(i), ',', num2str(j),')'));
                coeff = abs(corrcoef(exp(i,:),exp(j,:)));
                text(0.25*line(1), 0.9*line(2), strcat('X:',num2str(c1(1,2)),'   Y:',num2str(c2(1,2)),'   Co:',num2str(coeff(1,2))));
                print(strcat(figname, 'CCM', num2str(i),'to',num2str(j)),'-dpng');
            end
        end
    end
    [~, sort_order] = sort(all_values(:,1),'descend');
    sorted_all_values = all_values(sort_order,:);
    sorted_all_values = cat(2,sorted_all_values(:,2:end),sorted_all_values(:,1));
else
    if strcmp(technique,'mi') || strcmp(technique,'te')
        exp = cont_to_disc(exp,20);
    end
    for i=1:num_nodes
        for j=1:num_nodes
            if i==j
                continue;
            else
                i
                j
                if strcmp(technique,'co')
                    values = zeros(1, 2*max_lag + 1);
                    for k = 1:2*max_lag+1
                        s1 = align_series(exp(i, :), k - max_lag - 1, 21);
                        s2 = align_series(exp(j, :), max_lag + 1 - k, 21);
                        values(k) = abs(corr(s1', s2'));
                    end
%                     values = abs(xcorr(exp(i,:),exp(j,:),max_lag,'coeff'));
%                     values = values(:,max_lag+1:end);
                    for l=1:length(values)
                        if isnan(values(l))
                            values(l) = 0;
                        end
                    end
                elseif strcmp(technique,'gc')
                    values = granger_cause(exp(i,:),exp(j,:),1,max_lag);
                    if isnan(values)
                        values = 0;
                    end
                else
                    values = zeros(max_lag,1);
                    for k=1:max_lag
                        if strcmp(technique,'mi')
                            values(k) = mutual_information(exp(i,:),exp(j,:),k);
                        elseif strcmp(technique, 'te')
                            values(k) = transfer_entropy(exp(i,:),exp(j,:),k);
                        end
                    end
                end
                r = find(values==max(values));
                %all_values = cat(1,all_values,[values(1),i,j,r(1)]); %uncomment for max_lag test on CO
                all_values = cat(1,all_values,[max(values),i,j,r(1)]);
            end
        end
    end
    [~, sort_order] = sort(all_values(:,1),'descend');
    sorted_all_values = all_values(sort_order,:);
    sorted_all_values = cat(2,sorted_all_values(:,2:end),sorted_all_values(:,1));
end
[precision, recall, fscore, fpr, tpr] = fscore_roc_evaluation(num_nodes, grn, sorted_all_values);
end