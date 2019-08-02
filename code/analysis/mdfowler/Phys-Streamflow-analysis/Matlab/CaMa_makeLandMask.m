% This program creates a land mask for use with CaMa data - 1's are placed
% in grid cells that cover land (where there is data, at least), and 0's 
% are placed over the ocean. 
%
% Meg D Fowler, 2016-10-07

function landmask = CaMa_makeLandMask()

sp_file = ['/Volumes/MyPassport/CaMa/output/SP_cam3.5/outflwSP_cam3.5.nc'];

lon  = ncread(sp_file, 'lon');
lat  = ncread(sp_file,'lat');
sp_mean = ncread(sp_file,'ann_mean');%(m^3/s)

landmask = ones(size(sp_mean(:,:,1)));

ocean = find(isnan(sp_mean(:,:,1))==1);

landmask(ocean)=0;

% figure;
% contourf(lon,lat,landmask');      %It works :) 


end
