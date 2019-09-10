% extract morphological and intensity-based features from segmented 
% objects located within segmented cells
% 

function [spot_features, conn_feat_names, all_thumbnails, all_outlines, all_centroids] = ...
    object_feature_extraction( binary_map, intensity_map )

% defaults
thumbnail_size = 25; % odd positive integer, n x n thumbnail centered at weighted centroid
haralick_feats = [1 2 4:11 13]; % select features [1 14]
gauss_corr = 1; % always a single feature
haus_dist = 1; % always a single feature
conn_feat_names = parseNames({'Haralick','GaussCorr','HausDist','Morphological', 'Intensity', 'GradientIntensity'},haralick_feats,gauss_corr,haus_dist);
% ^ connects feature columns to names

% find connected components
cc = bwconncomp( binary_map );
obj = regionprops(cc, intensity_map, 'all'); 

% get centroids
all_centroids = cat(1,obj(:).WeightedCentroid) ;

% loop through each object
padding = 1 ; % for thumbnail image
spot_features = zeros( length(obj), length(conn_feat_names), 'single' ); 
all_thumbnails = cell( 1, length(obj) );
all_outlines = cell( 1, length(obj) );
for ii=1:length(obj)
    % make a thumbnail image
    bbox = obj(ii).BoundingBox ;
    xylims = [size(intensity_map,1), size(intensity_map,2)];
    [x1,x2,y1,y2] = get_bounding_box_corners(bbox, padding, xylims);
    
    thumbnail = intensity_map( y1:y2, x1:x2 ); % extracted spot with padding=1pix to help with haralick features
    
    % Haralick features
    glcm = graycomatrix( thumbnail, 'offset', [0 1], 'Symmetric', true);
    haralick = haralickTextureFeatures(glcm, haralick_feats);
    haralick = haralick( haralick_feats ) ;
                      
    % Distance to Gaussian spot   
    xx = 1:(x2-x1)+1; yy = 1:(y2-y1)+1;
    [X,Y] = meshgrid(xx,yy);
    Gsig = obj(ii).EquivDiameter/4 ; % radius ~ 2sigma
    G = 1/sqrt(4*pi*Gsig^2) .* exp( -( (X-(x2-x1+1)/2).^2+(Y-(y2-y1+1)/2).^2)./(4*Gsig^2) ) ;
    G = sum(obj(ii).PixelValues) .* G;  % scale intensities 
    Gcorr = corr2( intensity_map(y1:y2, x1:x2), G ); % correlation between images (one number)

    % Hausdorff distance
    bwG = zeros( size(G) ); 
    bwG( round( (size(bwG,1)+1)/2), round( (size(bwG,2)+1)/2 ) ) = 1; % centerpoint
    bwG = imdilate(bwG, strel('disk',round(obj(ii).EquivDiameter/2)));
    Hdist = calc_hausdorff( bwperim(bwG), bwperim(binary_map(y1:y2, x1:x2)), 99 ); % 99th percentile of distances for outlier robustness

    % store extracted features
    spot_features(ii,:) = [haralick', Gcorr, Hdist] ; 
    
    % store thumbnails for viewing later 
    o_center = all_centroids(ii,:);
    o_x = round(all_centroids(ii,1)-(thumbnail_size-1)/2:all_centroids(ii,1)+(thumbnail_size-1)/2) ;
    o_y = round(all_centroids(ii,2)-(thumbnail_size-1)/2:all_centroids(ii,2)+(thumbnail_size-1)/2) ;
    % ensure positive values
    while any(o_x<=0), 
        o_x=o_x+1; 
    end
    while any(o_y<=0),
        o_y=o_y+1;
    end
    % ensure within boundaries of image
    while any(o_x>=size(binary_map,1))
        o_x=o_x-1;
    end
    while any(o_y>=size(binary_map,1)) % assuming square image
        o_y=o_y-1;
    end
    all_thumbnails{ii} = intensity_map( o_y, o_x );
    all_outlines{ii} = bwperim( binary_map(o_y,o_x) );
    
    
%     [xx1,xx2,yy1,yy2] = get_bounding_box_corners(bbox, 5, xylims);
%     bigger_thumb = intensity_map( yy1:yy2, xx1:xx2 );
%     all_thumbnails{ii} = bigger_thumb;
%     bigger_overlay = binary_map( yy1:yy2, xx1:xx2 );
%     all_outlines{ii} = bwperim( bigger_overlay );
end
% Morphological features
morph_feats = [ [obj(:).Area]', ...
                [obj(:).Eccentricity]', ...
                [obj(:).EquivDiameter]', ...
                [obj(:).MajorAxisLength]', ...
                [obj(:).MinorAxisLength]', ...
                [obj(:).Orientation]', ...
                [obj(:).Perimeter]', ...
                [obj(:).Solidity]' , ...
                [obj(:).Perimeter]'.^2 ./ (4*pi*[obj(:).Area]'), ... % compactness
                [obj(:).Perimeter]' ./ cellfun(@(x) sum(sum(bwperim(x))), {obj(:).ConvexImage}, 'UniformOutput', true)', ... % roughness p/p_convex
                [obj(:).MajorAxisLength]' ./ [obj(:).MinorAxisLength]' ... % shape feature
                ];
morph_names = {'Area','Eccentricity','EquivDiameter',...
               'MajorAxisLength','MinorAxisLength','Orientation','Perimeter',...
               'Solidity','Compactness','Roughness','MajorToMinorAxisLengthRatio' }; 
conn_feat_names = [conn_feat_names morph_names(:)'];

% Intensity features
inten_feats = [ [obj(:).MaxIntensity]', ...
                [obj(:).MinIntensity]', ...
                [obj(:).MeanIntensity]', ...
                cellfun(@(x) var(x(:)), {obj(:).PixelValues}, 'UniformOutput', true)' ... % variance
               ];
inten_names = {'MaxIntensity','MinIntensity','MeanIntensity','VarianceIntensity'};
conn_feat_names = [conn_feat_names inten_names(:)'];

% Gradient intensity features
grad_obj = regionprops(cc, imgradient(intensity_map,'sobel'), 'all');
gradinten_feats = [ [grad_obj(:).MaxIntensity]', ...
                    [grad_obj(:).MinIntensity]', ...
                    [grad_obj(:).MeanIntensity]', ...
                    cellfun(@(x) var(x(:)), {grad_obj(:).PixelValues}, 'UniformOutput', true)' ... % variance
                   ];
gradinten_names = {'MaxGradIntensity','MinGradIntensity','MeanGradIntensity','VarianceGradIntensity'};
conn_feat_names = [conn_feat_names gradinten_names(:)'];


% combine features
spot_features = [spot_features, morph_feats, inten_feats, gradinten_feats];

end

% helper functions
function [hdist] = calc_hausdorff(ground_truth,test_mask,varargin)
% to be symmetric, both the forward and reverse are calculated
GT_edges = bwperim(ground_truth);
SEG_edges = bwperim(test_mask);
% Distance transform (default=Euclidean)
df = bwdist(GT_edges);
dr = bwdist(SEG_edges);
% edge distances only
ef = df( test_mask );
er = dr( ground_truth );
% symmetric calculation (forward and reverse)
if ~isempty(varargin)
    if varargin{1}<0 || varargin{1}>100
        error('OutOfBounds: Hausdorff distance percentile must be [0 100].');
    end
    k = varargin{1} ; % percentile for calculation
    hf = prctile(ef,k);
    hr = prctile(er,k);
else % standard to use max
    hf = max(ef);
    hr = max(er);
end
hdist = max( hf, hr );
end

function [x1,x2,y1,y2] = get_bounding_box_corners(bound_box, pad, lims)
    x1 = max( floor(bound_box(1))-pad,1) ;
    x2 = min( floor(bound_box(1))+ceil(bound_box(3))+pad, lims(1) );
    y1 = max( floor(bound_box(2))-pad,1) ;
    y2 = min( floor(bound_box(2))+ceil(bound_box(4))+pad, lims(2) );
end

function [feat_name_list] = parseNames( feat_names, haralick, gc, hd )

feat_name_list = [];
for ii=1:length(feat_names)
    switch feat_names{ii}
        case 'Haralick'
            for har=haralick
                feat_name_list = [feat_name_list, {strcat(feat_names{ii},num2str(har))} ];
            end
        case 'GaussCorr'
            feat_name_list = [feat_name_list, repmat(feat_names(ii),1,length(gc))] ;
        case 'HausDist'
            feat_name_list = [feat_name_list, repmat(feat_names(ii),1,length(hd))] ;
        case {'Morphological','Intensity','GradientIntensity'}
            % these are added later
        otherwise
            error( 'During extract_feature_from_objects, feature_extraction: Input feature name unknown.' );
    end
end

end

