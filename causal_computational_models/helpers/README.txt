Includes some helper functions. Also, some functions which can help in collating results and doing some post-evaluation analysis. See scripts:
output_aurocs.m: to tabulate AUROCs of all 5 DREAM4 networks of given size and for all PW metrics
plot_aggregate.m: to plot ROC curves of all 5 DREAM4 networks of given size and for all PW metrics on the same figure
plot_score_pdf.m: to plot the score distributions for all 5 DREAM4 networks of given size and for all PW metrics
anova_test.m: to do an ANOVA test of the score distributions for true and false edges for all 5 DREAM4 networks of given size
topk_edge_analysis.m: to find the node-degree distrbutions of topk links detected by the pairwise metrics