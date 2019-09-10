% Run principle components analysis
% 
% Copyright J. Hartman 2019, Cornell University


function [output_values,model_struct] = run_pca(data, var_thresh, options)

% Apply PCA to reduce correlation of object features and reduce dimension
[z_data,mu,stdev] = zscore(data); % standardize data parameters for future use
[coeffs,~,latent,~,explained] = pca( data, 'VariableWeights', 'Variance' );

% Remove noise using threshold and then tossing out some components.
% Here, we must be careful using PCA to actually reduce dimension before 
% clustering because there is no guarantee that variance of a feature is 
% equivalent to discrimination for classification.
pca_idx = cumsum(explained)<var_thresh ; % percent of variance retained
eigval = diag(latent(pca_idx));
eigvec = diag(stdev) \ coeffs(:,pca_idx); % PCA coefficients have been scaled 

% Optional sphering or whitening of the eigenvalues
if ~exist('options','var') % no input
    output_values = z_data * eigvec ; % classic PCA scores (not whitened)
elseif ~strcmp(options,'whiten') % input is not 'whiten'
    output_values = z_data * eigvec ; % classic PCA scores (not whitened)
else % PCA-whitening
    epsilon = 1e-5; % changing this value will smooth the output
    output_values = data * diag(1./sqrt(diag(eigval)+epsilon)) * eigvec ;
end

% Storing components for future use
model_struct.mu = mu;
model_struct.std = stdev;
model_struct.eigenvectors = eigvec;
model_struct.eigenvalues = eigval;
model_struct.explained = explained;

end
