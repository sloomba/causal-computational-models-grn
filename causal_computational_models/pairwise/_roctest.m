function [ precision, recall, fscore, fpr, tpr, sorted_all_values, all_values ] = roctest( exp, grn, technique, figname, max_lag )
%THIS SCRIPT WAS DEVELOPED IN THE EARLIER STAGES OF THE PROJECT. CAN BE IGNORED.
%Does ROC analysis using given pairwise technique, assuming
%num_nodes > num_time_steps.
%Finds optimal lag using AIC/BIC
%Doesn't count for self regulation
%techniques = 'ra', 'co', 'gc', 'mi', 'te', 'cm'
if nargin<5
    max_lag = 10;
    if nargin<4
        figname = '';
        if nargin<3
            technique = 'ra';
        end
    end
end
technique
warning off;
if (size(exp,1)>=size(exp,2))
    exp = exp';
end
num_nodes = size(exp,1);
% all_values = zeros(num_nodes*(num_nodes-1), 4);
all_values = [];
if strcmp(technique, 'mi') || strcmp(technique, 'te')
    exp = cont_to_disc(exp, 20);
end
if strcmp(technique,'ra')
    for i=1:num_nodes
        for j=1:num_nodes
            if i==j
                continue;
            else
                all_values(num_nodes*i + j, :) = [i, j, 0, 0];
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
            [~, ~ , ~ , SugiY , SugiX , origY , origX] = convergent_cross_map(exp(i,:)',exp(j,:)',1,4);
            SugiX(~isfinite(SugiX)) = 0;
            SugiY(~isfinite(SugiY)) = 0;
            c1 = corrcoef(SugiX, origX);
            c2 = corrcoef(SugiY, origY);
            all_values(num_nodes*(i-1) + (j-1), :) =  [c1(1,2), i, j, 1];	
            all_values(num_nodes*(j-1) + (i-1), :) = [c2(1,2), j, i, 1];
            %all_values = cat(1,all_values,[c1(1,2),i,j,1]);
            %all_values = cat(1,all_values,[c2(1,2),j,i,1]);
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
    for i=1:num_nodes
        for j=1:num_nodes
            if i==j
                continue;
            else
                i
                j
                if strcmp(technique,'co')
                    values = abs(xcorr(exp(i,:),exp(j,:),max_lag,'coeff'));
                    values = values(:,max_lag+1:end);
                    for l=1:length(values)
                        if isnan(values(l))
                            values(l) = 0;
                        end
                    end
                elseif strcmp(technique,'gc')
                    values = granger_cause(exp(i,:),exp(j,:),1,max_lag);
                    if isnan(values)
                        values = -Inf;
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
                % all_values(num_nodes*(i-1) + (j-1), :) = [max(values),i,j,r(1)]; % cat(1,all_values,[max(values),i,j,r(1)]);
                all_values = cat(1,all_values,[max(values),i,j,r(1)]);
            end
        end
    end
    [~, sort_order] = sort(all_values(:,1),'descend');
    sorted_all_values = all_values(sort_order,:);
    sorted_all_values = cat(2,sorted_all_values(:,2:end),sorted_all_values(:,1));
end
num_nodes, size(grn), size(sorted_all_values)
[precision, recall, fscore, fpr, tpr] = fscore_roc_evaluation(num_nodes, grn, sorted_all_values);
end