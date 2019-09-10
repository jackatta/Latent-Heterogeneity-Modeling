
% load viz if not already in workspace
if ~exist(viz)
    load( fullfile('data','extracted_viz_info.mat') );
end

% evaluate the clipped spots (need viz variable loaded)
figure,
clip_id = find( clipped.origidx );
viz.thumbnails = cat(2,viz.thumbnails{:});
viz.outlines = cat(2,viz.outlines{:});
nviz = 25; % make divisable by 5
r_id = randperm(length(clip_id),nviz);
r_id = clip_id(r_id);
for jj=1:length(r_id)
    subplot(nviz/5,5,jj),
    imshow( imOverlay( mat2gray(viz.thumbnails{r_id(jj)}), viz.outlines{r_id(jj)}, 'r' ) );
end

clear r_id nviz clip_id