%% Bilateral filtering for images (2D)
% [output]   =   BilateralFilter(I, W, SIGMA)
% 
% Standard nonlinear filtering technique usually used for preserving edges.
% This code is for monochromatic images, and is partly copied from
% BFILTER2 by Douglas Lanman on the MATLAB file exchange.
% 
% I : input 2D image (normalized) in format of double precision
% W : half-size of Gaussian bilateral filter window 
% SIGMA : 1x2 vector of bilateral filter standard deviations where SIGMA(1)
% is for the spatial domain and SIGMA(2) is for the intensity domain
% 
% Filter response is the product of a Gaussian-weighted distance function
% (SIGMA(1)) and a Gaussian-weighted intensity function (SIGMA(2)).
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Modified from BilateralFilter on MATLAB file exchange to use the function
% im2col, which can really speed up operations occurring in local windows.
% Local windows must be square [MxM] in the current implementation.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function [output]   =   BilateralFilter_fast(I, W, SIGMA)

% Verify that the input image exists and is valid
if ~isfloat(I) || 1 ~= size(I,3) || min(I(:)) < 0 || max(I(:)) > 1
   error(['Input image A must be a double precision ',...
          'matrix of size NxMx1 on the closed ',...
          'interval [0,1].']);      
end

% Pre-compute Gaussian distance weights
[X,Y]   =   meshgrid(-W:W,-W:W);
domfilt =   exp(-(X.^2+Y.^2)/(2*SIGMA(1)^2));

% prepare array for filtering
padSize = [W W];
padded_image = padarray( I, padSize, 'symmetric', 'both' );
col_image = im2col(padded_image,2*[W W]+1,'sliding'); % local 2D window becomes 1D column vector

% local calculations
halfway = ceil(size(col_image,1) / 2) ;
col_out = zeros(1, size(col_image,2)) ;
ii = 1 ;
while ii<=size(col_image,2) % while loops can be faster
% for ii=1:size(col_image,2)
    % apply intensity weights
    intfilt = exp(-(col_image(:,ii)- col_image(halfway,ii)).^2/(2*SIGMA(2)^2)) ;
    
    % calculate bilateral filter response
    BLfilt = intfilt .* domfilt(:);
    col_out(ii) = sum( BLfilt .* col_image(:,ii) ) ./ sum(BLfilt) ;
    
    ii = ii + 1 ;
end

% reshape back into image
[MM,NN] = size(padded_image);
M = 2*W+1; N = 2*W+1;
output = reshape(col_out,MM-M+1,NN-N+1);


end