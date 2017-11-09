function [p_vals] = anova_test(num_nodes, pw_or_ige)
%Do an ANOVA test on score distributions. Example usage:
%anova_test(10,'pw');
    if nargin<2
        pw_or_ige = 'pw';
    end
    p_vals = [];
    for idx=1:5
        dataset_name = strcat('dream_', num2str(num_nodes),'_',num2str(idx),'_correct');
        load(strcat('../data/data_', dataset_name), 'grn'); % Change path of data file for grn here
        grn = grn(:, 1:2);
        load(strcat('../data/result_', dataset_name, '_', pw_or_ige));  % Change path for result file here
        names = {};
        count_t = length(grn);
        for j=1:count_t
            names{j} = 'true';
        end
        count_f = num_nodes*(num_nodes-1) - count_t;
        for j=count_t+1:count_t+count_f
            names{j} = 'false';
        end
        anova_results = {};
        curr_p = [];
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
                true_edges = values(index, 4)';
                index = find(~r);
                false_edges = values(index, 4)';
                groups = [true_edges, false_edges];
                [p, anova_table, stats] = anova1(groups, names);
                close(); close();
                anova_result{i}.p = p;
                curr_p = cat(1, curr_p, p);
                anova_result{i}.anova_table = anova_table;
                anova_result{i}.stats = stats;
            end
            p_vals = cat(2, p_vals, curr_p);
        else
            %% Code for IGE
            % TODO
        end
        save(['../data/anova_' dataset_name '_' pw_or_ige], 'anova_result');
    end
end