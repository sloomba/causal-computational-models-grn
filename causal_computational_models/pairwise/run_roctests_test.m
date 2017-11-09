function [ pw_out ] = run_roctests( exp, grn, figname, resname)
%RUN_ROCTESTS for time-lag PW analysis
method = containers.Map(1,'Random');
method_code = containers.Map(1,'ra');
method(2) = 'CO';
method_code(2) = 'co';
method(3) = 'GC';
method_code(3) = 'gc';
method(4) = 'MI';
method_code(4) = 'mi';
method(5) = 'TE';
method_code(5) = 'te';
method(6) = 'CCM';
method_code(6) = 'cm';
test_for = [2,3,6]
for testid = 1:length(test_for)
    j = test_for(testid);
    precision = {[] [] [] [] [] [] [] [] [] [] [] [] [] [] [] []};
    recall = {[] [] [] [] [] [] [] [] [] [] [] [] [] [] [] []};
    fscore = {[] [] [] [] [] [] [] [] [] [] [] [] [] [] [] []};
    fpr = {[] [] [] [] [] [] [] [] [] [] [] [] [] [] [] []};
    tpr = {[] [] [] [] [] [] [] [] [] [] [] [] [] [] [] []};
    values = {[] [] [] [] [] [] [] [] [] [] [] [] [] [] [] []};
    unsorted_values = {[] [] [] [] [] [] [] [] [] [] [] [] [] [] [] []};
    if nargin>=3
        for i=1:16
            [precision{i}, recall{i}, fscore{i}, fpr{i}, tpr{i}, values{i}, unsorted_values{i}] = roctest(exp,grn, method_code(j), i);
            figure(2); clf;
            subplot(1,2,1);
            hold;
            plot(precision{i}, 'color', 'r');
            plot(recall{i}, 'color', 'g');
            plot(fscore{i});
            legend('Precision', 'Recall', 'f-score');
            xlabel('Links');
            ylabel('Score');
            hold;
            title(strcat('Precision, Recall and f-score plot for:', method(j),' MAX_LAG=',num2str(i)));
            subplot(1,2,2);
            plot(fpr{i}, tpr{i});
            xlabel('False Positive Rate');
            ylabel('True Positive Rate');
            text(0.5, 0.2, num2str(abs(trapz(fpr{i}, tpr{i}))));
            title('ROC plot');
            print(strcat(figname, method(j), num2str(i), 'ROCplotPW'),'-dpng');
            figure(2); clf;
            hold;
            scatter(grn(:,1),grn(:,2),'bo');
            temp = values{i};
            scatter(temp(1:size(grn,1),1),temp(1:size(grn,1),2),'g*');
            legend('Ground', 'Found');
            xlabel('From node');
            ylabel('To node');
            hold;
            title(strcat('Scatter plots for causal relations:', method(j),' MAX_LAG=',num2str(i)));
            print(strcat(figname, method(j), num2str(i), 'ScatterPW'),'-dpng');
        end
    end
    pw_out = {precision, recall, fscore, fpr, tpr, values, unsorted_values};
    save(strcat(resname,'lagtest',method(j)),'pw_out');
end