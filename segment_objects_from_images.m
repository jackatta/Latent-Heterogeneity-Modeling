% Segment endosome-like objects from images using Spot Sauvola

function [segmented_masks,original_images,spots] = segment_objects_from_images(img_labels, cell_masks, img_dir, seg_dir)

cell_areas = squeeze(sum(sum(cell_masks,1),2));

% initialize some variables
tnim = height(img_labels); % total num images
segmented_masks = false(2048,2048,tnim);
original_images = zeros(2048,2048,tnim, 'uint16');
spots.data = table();

% go through and process the images
nim = 1; % total image counter
tif_list = dir( fullfile(img_dir,'*.tif') ); % NOTE that tiff files were ordered manually in the folder
tic
for ii=1:length(tif_list)
    im_info = imfinfo( fullfile(img_dir,tif_list(ii).name) );
    for jj=1:length(im_info)
        % load image
        original_images(:,:,nim) = imread( fullfile(img_dir,tif_list(ii).name), jj, 'Info', im_info );
        
        % segment objects 
        Spot_Sauvola;
        mask_name = strrep(tif_list(ii).name,'.tif','_classified_proposed.tif') ;
        imwrite( final_image, fullfile(seg_dir,mask_name), 'WriteMode', 'append', 'Compression', 'none' );
        segmented_masks(:,:,nim) = final_image;
        
        % extract info from objects
        [temp_features, temp_names, viz.thumbnails{nim}, viz.outlines{nim}, ~] = ...
          object_feature_extraction( final_image, RescaleImage( original_images(:,:,nim), 2^16-1 )) ;
        
        % attach image info to objects
        temp_features = abs(temp_features); % Haralick13 complex number in 2 cases... not sure why
        temp_table = array2table(temp_features,'VariableNames',temp_names);
        temp_table.CellArea = repmat(cell_areas(nim),height(temp_table),1);
        for kk=2:length(img_labels.Properties.VariableNames) % don't bother with FileName
            temp_table.(img_labels.Properties.VariableNames{kk}) = repmat(img_labels{nim,kk},height(temp_table),1);
        end
        temp_table.UniqIDperImage = (1:height(temp_table))';
        spots.data = vertcat(spots.data,temp_table);
        
        if rem(nim,10)==0
            fprintf('%3.0f/%3.0f complete, %3.1fsec',nim,tnim,toc)
        end
        nim=nim+1;
    end
end
spots.data.UniqID = (1:height(spots.data))';

% also save as mat file for quicker loading if needed
save( fullfile('data','extracted_viz_info_09092019.mat'), 'viz', '-v7.3' );

end
