function [] = plot_score_pdf(num_nodes, pw_or_ige)
%Find score distributions. Example usage: plot_score_pdf(10,'pw');
    if nargin<2
        pw_or_ige = 'pw';
    end
    for idx=1:5
        dataset_name = strcat('dream_', num2str(num_nodes),'_',num2str(idx),'_correct');
        load(strcat('../data/data_', dataset_name), 'grn'); % Change path of data file for grn here
        grn = grn(:, 1:2);
        load(strcat('../data/result_', dataset_name, '_', pw_or_ige));  % Change path for result file here
        if(strcmp('pw', pw_or_ige)) % For pw results
            %% Code for PW
            scores = pw_out{6};
            method = containers.Map(1,'Random');
            method(2) = 'CO';
            method(3) = 'GC';
            method(4) = 'MI';
            method(5) = 'TE';
            method(6) = 'CCM';
            for i=2:6
                values = scores{i};
                [~, r] = ismember(values(:, 1:2), grn, 'rows');
                index = find(r);
                true_edges = values(index, 4);
                [y_all, x_all] = ksdensity(values(:, 4));
                [y_tru, x_tru] = ksdensity(true_edges);
                plot(x_all, y_all);
                hold all;
                plot(x_tru, y_tru);
                legend('all edges','true edges')
                xlabel('Score');
                ylabel('pdf');
                title(['Probability Distribution of edge scores for ' method(i)]);
                print(['../plots/plots_analysis' dataset_name '_' method(i) 'pdf'], '-dpng');
                close();
            end
        else
            %% Code for IGE
            % TODO
        end
    end
end