% Hypothesis testing on macropinocytosis activity scores
% Inference performed using Wald tests on nonlinear trends using the input
% weights as the trend parameters.

function [p_table] = activity_score_inference(mdl,wgt,opt)
%% Determine trend parameters from input wgts
% mdl coefficients can be in different order than weights b/c reference level is removed from model
[emmM]      = emmeans(mdl,{opt.gmm.latent_var},'balanced'); % just getting order of mpc_class for model coefficients
char2num    = double(char(emmM.table.(opt.gmm.latent_var)));
row2mdl     = char2num-min(char2num)+1;

% rank and order the weights
[~,row2rank] = sort(wgt); % rank is ascending
row_wgt = wgt./min(wgt);
[~,~,rank2mdl] = intersect( row2mdl, row2rank, 'stable');
mdl_wgt     = row_wgt(row2rank(rank2mdl));
rank_class  = cellstr(char('A'-1+row2rank));

m_trends = mdl_wgt ./ sum(mdl_wgt); % make it a probability

%% Null statistical hypothesis tests
p_table = trend_hypothesis_testing(mdl,m_trends,opt);

end
