% Extract objects using proposed method based on the marker-controlled
% watershed transformation, where the markers are obtained from a sensitive
% segmentation (h-dome) reconstructed under a coarser segmentation.
% 

starting_image = double(original_images(:,:,nim)) ;
      
rad = 2 ;
smooth_image = BilateralFilter_fast( RescaleImage(starting_image,2^16), rad, [2 2] );

h=0.1;
rad_h = 5 ; % of flat disk structuring element
[enhanced_image,bkgd] = HDomeTransform( smooth_image, h, rad_h, 'recon conn' );

% low FN but low FP and cannot conserve area 
sauv_win = [25 25];
sauv_thresh = 0.005 ;
enhanced_mask = sauvola_inverse( enhanced_image, sauv_win, sauv_thresh) ; % from the FEX and nearly matches SauvolaThreshold written manually just wayyyy faster

% high FP and conserve area but low FN
unenhanced_mask = sauvola_inverse( smooth_image, sauv_win, sauv_thresh );

% combining locally thresholded masks to reduce FP
combined_masks = imreconstruct( enhanced_mask, unenhanced_mask );
combined_masks = imfill( combined_masks, 'holes' );

% further reduce FP with morphological cleanup
dom.shape = 'disk'; 
dom.size = 1 ;
spot_markers = imerode(combined_masks, strel(dom.shape,dom.size));

% watershed for accurate spot areas
[flooded_image,~] = MarkerSeededWatershed(spot_markers,smooth_image);

% add a quick classifier on top
% include only objects within the identified cell boundaries
binary_map = flooded_image & cell_masks(:,:,nim);
% extract a couple features from connected components
cc = bwconncomp(binary_map);

if cc.NumObjects~=0
    grad_stats = regionprops( cc, imgradient(smooth_image,'sobel'), 'MeanIntensity' ); 
    avg_grad = [grad_stats(:).MeanIntensity];
    int_stats = regionprops( cc, smooth_image, 'MeanIntensity' );
    avg_int = [int_stats(:).MeanIntensity] ;
    % cluster
    ktest = [avg_grad' avg_int'];
    nclust = min(cc.NumObjects,20); % foreground and background only
    [vec_LM,centers] = kmeans(ktest, nclust, 'Replicates', 10 ); % label matrix result from kmeans algorithm
    background_id = find( centers(:,1)==min(centers(:,1)) ) ;
    vec_LM = vec_LM~=background_id ;

    % update region labels
    final_image = ismember( labelmatrix(cc), find(vec_LM==1) ) ;
    n_obj = sum(vec_LM==1);
else
    final_image = false(size(binary_map));
    n_obj = 0;
end

% figure, imshowpair( smooth_image, final_image )
% 
% 
% overlay_image = imOverlay( adapthisteq(RescaleImage(starting_image)), bwperim(spot_markers), 'r', bwperim(flooded_image), 'b', bwperim(final_image), 'g' );
% figure, imshow( overlay_image )
% figure, imshow( adapthisteq(RescaleImage(starting_image)), [])
