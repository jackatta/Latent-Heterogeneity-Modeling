% Null hypothesis testing of ordinal trends. Trends determined by a ranked
% ordering latent object classes where the ranking is based on a single
% feature value.
% Inference performed using Wald tests.

function [p_table,mpc_class_order] = object_feature_inference(mdl,objdata,trendname,opt)

%% Determine ordering for single feature

% get robust ordering metric for endosome classes
mpc_response = objdata(:,{trendname, opt.gmm.latent_var});
mpc_response = grpstats(mpc_response,{opt.gmm.latent_var},{'median'},'Datavars',trendname); % stats of response var per class

% set up ordering of mpc classes for trend testing
[~,mpc_response.ind_order] = sort( mpc_response.(sprintf('median_%s',trendname)) );
orderstat(mpc_response.ind_order,1) = 1:height(mpc_response);
orderspace = orderstat ./ min(orderstat); % spaced mpc class trend
mpc_response.order_trend = orderspace - mean(orderspace); % centered

% mdl coefficients can be in different order than weights b/c reference level is removed from model
[emmM] = emmeans(mdl,{opt.gmm.latent_var},'balanced'); % just getting order of mpc_class for model coefficients
char2num = double(char(emmM.table.(opt.gmm.latent_var)));
mpc_response.mdlcoefs = char2num-min(char2num)+1;
[~,~,mpc_response.mdl_order] = intersect( mpc_response.mdlcoefs, mpc_response.ind_order, 'stable');
mpc_response.mdl_trend = mpc_response.order_trend(mpc_response.ind_order(mpc_response.mdl_order));
mpc_class_order = sortrows([emmM.table.(opt.gmm.latent_var) num2cell(mpc_response.mdl_order)],2);
mpc_class_order = mpc_class_order(:,1)';

%% Null statistical hypothesis testing
p_table = trend_hypothesis_testing(mdl, mpc_response.mdl_trend, opt);

end