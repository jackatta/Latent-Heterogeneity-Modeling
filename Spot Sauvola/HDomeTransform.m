% The  h-dome transformation performs a geodesic dilation until stability 
% (i.e. opening by reconstruction) where I-h is the marker and 
% I is the mask.
% 
% By J. Hartman 6/2017
% 
% This Matlab function is based off of part of the code package released as a
% supplementary material for the article: 
%"Evaluation of methods for detection of fluorescence labeled subcellular
% objects in microscope images" by Ruusuvuori
% Website: http://www.cs.tut.fi/sgn/csb/subcell

function [hdImage,bkgd] = HDomeTransform(I, h, discradius, option)

switch option
    case 'recon conn'
        bkgd      =   imhmax( I, h ) ; % use default 8-point connectivity
    case 'recon se'
        % the purpose of this option is if user wants to use se with
        % dimensions greater than 3x3. For example, an 11x11 disk.
        se      =   strel( 'disk', discradius ) ;
        marker  =   I-h ;
        bkgd      =   ReconstructionWithSE( marker, I, se ) ;
end

hdImage     =   I - bkgd ;

end

% Because this function doesn't use MATLAB's native function, and it is 
% therefore not written according to method in Vincent 1993 nor has it been
% as MEX-file, it runs MUCH slower. 

function [output]   =   ReconstructionWithSE( marker, mask, se )

J       =   double(marker) ;
mask    =   double(mask) ;
flag    =   0 ;
tempI   =   -Inf * zeros( size(J), 'double' ) ;
convergencecriteria     =   1e-7 ;
disp('Running grayscale reconstruction with structuring element');
counter = 1 ;
while flag==0
	tempJ       =   min( imdilate(J, se), mask ) ;
	residual(counter)    =   sum( sum(abs(J-tempJ)) ) ;
	if residual(counter) < convergencecriteria
        flag = 1 ;
    end
            
	tempI       =   max( tempI, tempJ ) ;
	J           =   tempJ ; % update for next iteration
    counter = counter + 1 ;
%     if counter==50
%         flag = 1 ;
%     end
end
  
output  =   tempI ;

end
