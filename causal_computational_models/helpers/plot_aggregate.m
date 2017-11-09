function out = plot_aggregate(node)
%Aggregate ROC plots. Example usage: plot_aggregate(10);
prefix = '../data/result_dream';
suffix = 'correct_pw'; %only PW analysis right now
    data = {{} {} {} {} {}};
    for i=1:5
        data{i} = load(strcat(prefix, '_', num2str(node), '_' ,num2str(i), '_', suffix));
        data{i} = data{i}.pw_out;
    end
    data{6} = {};
    init = zeros(node*(node-1), 1);
    for i=1:5
        data{6}{i} = {init init init init init init};
    end 
    for i=1:5
        for j=1:5
            for k=1:6
            data{6}{j}{k} = data{6}{j}{k} + data{i}{j}{k}/5;
            end
        end
    end
    out = data{6};
    save(strcat(prefix, '_aggregated_', num2str(node), '_', suffix),'out');
    clear data;
    %Precision Recall Curve plot
    hold;
    plot(out{2}{2}, out{1}{2}, 'blue');
    plot(out{2}{3}, out{1}{3}, 'g');
    plot(out{2}{4}, out{1}{4}, 'r');
    plot(out{2}{5}, out{1}{5}, 'y');
    plot(out{2}{6}, out{1}{6}, 'black');
    hold;
    ylabel('Precision');
    xlabel('Recall');
    title(strcat('Precision Recall curve average over 5 DREAM4 networks of size ', num2str(node)));
    legend('Correlation', 'Granger Causality', 'Mutual Information', 'Transfer Entropy', 'Convergent Cross Map');
    print(strcat('../plots/AggPRC_dream_', num2str(node), suffix), '-dpng');
    close();
    %ROC curve plot
    hold;
    plot(out{4}{2}, out{5}{2}, 'blue');
    plot(out{4}{3}, out{5}{3}, 'g');
    plot(out{4}{4}, out{5}{4}, 'r');
    plot(out{4}{5}, out{5}{5}, 'y');
    plot(out{4}{6}, out{5}{6}, 'black');
    hold;
    ylabel('True Positive Rate');
    xlabel('False Positive Rate');
    title(strcat('ROC curve average over 5 DREAM4 networks of size ', num2str(node)));
    legend('Correlation', 'Granger Causality', 'Mutual Information', 'Transfer Entropy', 'Convergent Cross Map', 'Location', 'SouthEast');
    print(strcat('../plots/plots_analysis/AggROC_dream_', num2str(node), suffix), '-dpng');
    close();
end