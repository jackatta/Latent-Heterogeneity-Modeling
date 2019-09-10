
clear clipped

objectdata = spots.data(:,1:32);

% clip data too high using quantile threshold
clipped.names = {'HausDist','Area','Eccentricity','EquivDiameter','MajorAxisLength','MinorAxisLength','MajorToMinorAxisLengthRatio'};
clipped.qntile_thresh = 0.9995 ; % as a percent for quantile (cdf)
clipped.origidx = false(1,height(objectdata));
for ii=length(clipped.names):-1:1
    [c,e] = histcounts( objectdata{:,clipped.names{ii}}, 5000 );
    h_cdf = cumsum(c);
    h_cdf = h_cdf ./ max(h_cdf);
    clipped.origidx( objectdata{:,clipped.names{ii}} > min(e(h_cdf>clipped.qntile_thresh)) ) = true;
%     figure, plot(h_cdf);
%     title(clipped.names{ii});
%     figure, histogram(objectdata{:,clipped.names{ii}}, 5000);
%     title(clipped.names{ii});
end

% clip spots too small or anti-correlated to Gaussian spot
clipped.origidx( objectdata{ :, 'Area' } < 3 ) = true  ;
clipped.origidx( objectdata{ :, 'GaussCorr' } < -0.1 ) = true;

fprintf('Removed %3.0d outlier data points out of %3.0d total objects during clipping \n', sum(clipped.origidx), height(objectdata) );

clipped.data = spots.data(clipped.origidx,:);
spots.clipped = spots.data(~clipped.origidx,:);

clear objectdata c e h_cdf
