% Weighted marginalization of latent endosome class

function [m_emm, prob_mpc] = get_pheno_score(mdl, objdata, opt)
%% Define weights for probability of true macropinosome (big & bright)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % This section would vary for different biological applications
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Size
wgt_var{1} = 'EquivDiameter'; % this must be object feature name
pix2micron = 100/966; % um/pix as measured from scale with microscope settings
gam_shape = 6; gam_scale = 1000/gam_shape; % subjectively determined assuming 200nm+ with long tail
k{1} = makedist('Gamma','a',gam_shape,'b',gam_scale); % probability density of mpc given diameter
% figure, plot( (10:10:2000), pdf(k{1},(10:10:2000)) );

% Intensity
wgt_var{2} = 'MeanIntensity'; % this must be object feature name
k{2} = makedist('Logistic','mu',0.2,'sigma',0.04); % the cdf of k{2} is the probability density of mpc given intensity
% figure, plot( (0:0.01:1), cdf(k{2},(0:0.01:1)) );

% Total probability for each class
prob_perclass = grpstats(objdata(:,[wgt_var(:)',opt.gmm.latent_var]),{opt.gmm.latent_var},{'median'},'Datavars',[wgt_var(:)]'); % stats of response var per mpc_class
wgt = pdf(k{1},prob_perclass.(sprintf('median_%s',wgt_var{1}))*pix2micron*1000) ...
   .* cdf(k{2},prob_perclass.(sprintf('median_%s',wgt_var{2})));
prob_mpc = wgt ./ sum(wgt); % normalize to be probability

%% Form output emm table for the knowledge-weighted scores

% Estimated marginalized mean (random effects marginalized but not latent endosome class)
emm = emmeans(mdl,'unbalanced','effects','asfit','Offset', opt.glme.emm_offset);

% Marginalize latent endosome class using weights
bio_grp = mdl.Formula.FELinearFormula.PredictorNames;
bio_grp(strcmp(bio_grp,opt.gmm.latent_var)) = []; % remove latent variable that is marginalized
m_emm.table = grpstats(objdata,bio_grp,{'numel'},'datavars','CellArea'); % just getting table structure, don't need this actual data
m_emm.table(:,end-1:end) = []; % remove GroupCount and numel_CellArea
m_emm.table{:,{'Estimated_Marginal_Mean','SE','CI'}} = zeros(height(m_emm.table),3); % initializing new columns
for ii=1:height(m_emm.table)
    idx = ismember(emm.table(:,bio_grp), m_emm.table(ii,bio_grp));
    m_emm.table.Estimated_Marginal_Mean(ii) = emm.table.Estimated_Marginal_Mean(idx)' * prob_mpc;
    m_emm.table.SE(ii) = emm.table.SE(idx)' * prob_mpc;
end

end
