% Latent heterogeneity modeling first trains an object classifer that uses
% learned semantic dimensions and then fits a model that uses known
% experimental design
% 
% datatable is class table with all feature data and category labels
% data_flag is logical vector to provide table columns that contain data
%   rather than labels of the data
% strat_labels is cell array of strings to provide table column names for
%   variables to be used in stratification
% knw_labels is cell array of strings to provide table column names for
%   variables to be used as fixed effects (in addition to latent var)

function [mdl, datatable, pca_scores] = LatentHeterogeneityModeling(datatable, data_flag, opt)

objdata = datatable{:,data_flag};

%% Unsupervised training of object classifier (PCA + GMM)
train_and_apply_classifier; % learns mdl and adds column (latent_name) to datatable
save(fullfile('data','object_classification_model.mat'), 'mdl');

%% (Supervised) Fitting of phenotype classifier (GLME)
classcount_regression; % fits mdl.glme
save(fullfile('data','object_class_count_model.mat'), 'mdl.glme');

end