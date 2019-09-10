strat_datatable = datatable(:,opt.pca.strat_vars);

%% Stratify training subset
% Use a stratified, random subset of the feature data to train the model in 
% order to improve efficiency without loss of accuracy
mdl.objtraining.frac = opt.pca.strat_frac; % percent of data for training
mdl.objtraining.stratgroup = unique( strat_datatable );
mdl.objtraining.idx = false(1,height(datatable));
for ii=1:height(mdl.objtraining.stratgroup)
    idx = find(ismember(strat_datatable,mdl.objtraining.stratgroup(ii,:))); % this is really slow, was considerably faster doing nester for loops
    nobs = length( idx ); 
    mdl.objtraining.idx( idx(randperm(nobs,ceil(nobs*mdl.objtraining.frac))) ) = true; 
end

mdl.objtraining.data = datatable( mdl.objtraining.idx, : );

%% PCA
[mdl.objtraining.pca_scores, mdl.pca] = run_pca( double(objdata(mdl.objtraining.idx,:)), opt.pca.var_thresh );

%% GMM
% Initialize GMM training using centroids from kmeans clustering 
cInd = kmeans( mdl.objtraining.pca_scores, opt.gmm.n_clusters, 'EmptyAction', 'singleton', 'Replicates', 10, 'MaxIter', 300 );
mdl.gmm = fitgmdist( mdl.objtraining.pca_scores, opt.gmm.n_clusters, ...
                   'Start', cInd, ... % 'plus' seeds clusters with kmeans++
                   'Replicates', 1, ... % if 'plus', set to 10+
                   'CovarianceType', 'full', ...
                   'SharedCovariance', false, ...
                   'RegularizationValue', 1e-5, ...
                   'Options', statset('MaxIter', 1000) );
mdl.objtraining.gmm_labels = cluster(mdl.gmm, mdl.objtraining.pca_scores); % apply labels to training data


%% Apply classifier to all objects

% Standardize data using same mean and stdev as training data
zdata = double(objdata);
zdata = bsxfun(@minus,zdata,mdl.pca.mu);
zdata = bsxfun(@rdivide,zdata,mdl.pca.std);

% Project to PCA
pca_scores = zdata * mdl.pca.eigenvectors;

% Classify with GMM
temp_label = cluster( mdl.gmm, pca_scores );
datatable.(opt.gmm.latent_var) = cellstr(char('A'-1+temp_label)); % make classes alphabet

clear strat_datatable cInd ii idx nobs zdata temp_label