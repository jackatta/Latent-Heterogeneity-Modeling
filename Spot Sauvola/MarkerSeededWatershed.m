% Marker-seeded watershed using the distance transform of the markers

function [bwOutput,ws_lines] = MarkerSeededWatershed( bw_markers, im )
% make distance map and skeleton by influence zone map from markers
m_dist = bwdist(bw_markers) ; %distance transform
skiz   = watershed(m_dist)==0 ; %skeleton by influence zones
skizd  = imdilate(skiz,ones(3,3)); % this ensures that seeds don't touch skiz
newmarkers = bw_markers & ~skizd;
 
% use image gradient to find edges
[gradmag,~] = imgradient( im, 'sobel' ); % gradient magnitude
gradmag2 = imimposemin( gradmag, newmarkers | skiz, 4 ) ;

% perform seeded watershed
bwOutput = watershed(gradmag2) ;
ws_lines = bwOutput==0;
bwOutput = imclearborder(bwOutput)>0;

%         overlay = imdilate(Edges,ones(2,2)) | seeds | skiz ;
%         figure, imshowpair(Rescale(toanalyze), overlay*0.3, 'Scaling', 'joint'  );

end