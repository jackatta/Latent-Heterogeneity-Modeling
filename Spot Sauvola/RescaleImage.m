% rescale image to floating point range [0 1]
% default is to use max == 1 (equal to mat2gray) 
% or can provide the desired normalization value

function [I]   = RescaleImage( I, varargin )

I   = double(I);
I   = I - min(I(:)) ;

if nargin<2
    plotopt = 'noplot';
    max_val = max(I(:)) ;
elseif ischar(varargin{1})
    plotopt     = varargin{1};
    strTitle    = varargin{2};
    max_val = max(I(:)) ;
else
    plotopt = 'noplot';
    max_val = varargin{1} - min(I(:));
end


I   = I ./ max_val ;

if strcmp(plotopt,'plot')
    figure, imshow(I,[]); title(strTitle); 
end

end