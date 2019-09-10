% Spot Sauvola and Latent Heterogeneity Modeling for macropinocytosis
% scoring and characterization
% 
% Copyright J. Hartman 2019, Cornell University
% 
% Citation: John Hartman, Andrea Ortuno Urquidi, Fredrik Thege, Brian J
% Kirby. "Spot Sauvola and Latent Heterogeneity Modeling for scoring and 
% profiling macropinocytosis". Cytometry Part A (2019)

clear

add_scripts_to_path

[img_labels] = define_image_order; % Manually define image order for stacks

% Upload cell masks segmented previously from in-focus transmitted light
% images using the Trainable Weka Segmentation plugin in ImageJ/FIJI
load( fullfile('segmentations','Transmitted light channel','cell_masks.mat') );


%% Segment spots using Spot Sauvola and extract their features
segmentspots = true; % true = run segmentation, false = load data

% tell where files would be and what to name them
img_dir = fullfile('images','TMR channel');
img_name = fullfile(img_dir,'original_images.mat');
seg_dir = fullfile('segmentations','TMR channel');
seg_name = fullfile(seg_dir,'spot_masks.mat');
obj_name = fullfile('data','spots.mat');
if segmentspots
    [spot_masks, original_images, spots] = segment_objects_from_images(img_labels, cell_masks, img_dir, seg_dir);
    save( seg_name, 'spot_masks', '-v7.3' );
    save( img_name, 'original_images', '-v7.3' );
    save( obj_name, 'spots', '-v7.3' );
    clip_objects;
else
    load( seg_name ); % load segmentation maps
    load( img_name ); % load raw images
    load( obj_name ); % load object table
    clip_objects;
end
% % (optional: inspect clipped objects)
% inspect_clipped;

%% Remove unnecessary data
% the DEX- treatment condition was a helpful control but does not need to be analyzed
spots.clipped( ismember(spots.clipped.Treatment,'DEX-'), : ) = [];
cell_masks(:,:, ismember(img_labels.Treatment,'DEX-')) = [];
img_labels( ismember(img_labels.Treatment,'DEX-'), : ) = [];

%% Latent heterogeneity modeling to profile macropinocytosis
% Define optional parameters
options.pca.strat_vars = {'CellLine','Experiment','Treatment'}; % variables whose combinations define stratification groups
options.pca.strat_frac = 0.2; % percent of objects sampled per stratification group
options.pca.var_thresh = 99.9; % percent variance retained
options.gmm.n_clusters = 10; % number latent classes of endosomes
options.gmm.latent_var = 'EndosomeClass'; % variable name for latent class
options.glme.img_vars = img_labels.Properties.VariableNames(2:end); % variables whose combinations define individual images
options.glme.fixed_effects = {'CellLine','Treatment'}; % variables that are fixed effects (not including latent class which is added automatically)
options.glme.random_effects = ' + (1|Experiment) + (1|Experiment:Treatment) + (1|Experiment:Treatment:TechRep)'; % random effects formula
options.glme.offset = squeeze( sum(sum(cell_masks,1),2) ); % offset for each image
options.glme.emm_offset = log(mean(options.glme.offset)); % re-introduced offset to obtain final estimates per image
% Do the modeling
[mdl_LHM, spots.clipped, spots.pca_scores] = ...
    LatentHeterogeneityModeling(...
        spots.clipped, ... % table
        [true(1,32) false(1,6)], ... % which cols are data
        options);

%% Calculate Macropinocytosis Activity Score using weighted marginalization
[Activity_Score, MPC_wgts] = get_pheno_score(mdl_LHM.glme, spots.clipped, options);


%% Hypothesis testing

% Activity Score
hypothesis_tables.scores = activity_score_inference(mdl_LHM.glme, MPC_wgts, options);

% Trends of features for discovery (mult hypothesis correction req'd)
test_features = spots.clipped.Properties.VariableNames(1:32);
for ii=1:length(test_features)
    newfield = [test_features{ii} 'trend'];
    [temp_table, trend_order.(test_features{ii})] = object_feature_inference(mdl_LHM.glme,spots.clipped,test_features{ii},options);
    % only want to get p-values for RhoG trends because specific to MPC
    hypothesis_tables.discovery.(newfield) = temp_table(:,'linear');
end


%% Plotting
% Plot Scores
Activity_Score = relevel_table(Activity_Score, {'CellLine','Treatment'}, {'BxPC3','BxGR080C','BxGR360C'}, {'EGF-','BASELINE','VEH','EIPA40','EIPA80'});
h.marg.bycell = emmip( Activity_Score, '~ Treatment | CellLine' );
h.marg.compcel = emmip( Activity_Score, 'CellLine ~ Treatment' );

% sample code but can change class ordering as desired
Explore_Subpop = emmeans(mdl_LHM.glme,'unbalanced','effects','asfit','Offset', log(mean(options.glme.offset)));
Explore_Subpop = relevel_table(Explore_Subpop, {'EndosomeClass','CellLine','Treatment'}, {'A','B','C','D','E','F','G','H','I','J'}, {'BxPC3','BxGR080C','BxGR360C'}, {'EGF-','BASELINE','VEH','EIPA40','EIPA80'});
h.bycell = emmip(Explore_Subpop, 'EndosomeClass ~Treatment | CellLine');
