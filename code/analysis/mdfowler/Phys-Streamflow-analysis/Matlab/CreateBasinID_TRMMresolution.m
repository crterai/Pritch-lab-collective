% Create and save an array with Basin ID values, in order to separate out
% basins in CaMa data. 
% 
% Key: 
%   1       Yukon
%   2       Mackenzie
%   3       Columbia
%   4       Nelson
%   5       Mississippi
%   6       St-Lawrence
%   7       Rhine
%   8       Danube
%   9       Dniepr
%   10      Volga
%   11      Ob
%   12      Yenisei
%   13      Lena
%   14      Amur
%   15      Orinoco
%   16      Parana
%   17      Amazonas
%   18      Niger
%   19      Congo
%   20      Nile
%   21      Zambezi
%   22      Tigris & Euphrates
%   23      Indus
%   24      Ganges & Brahmaputra
%   25      Mekong
%   26      Huang-he
%   27      Yangtze
%   28      Murray & Darling
%   29      Don
%   30      Fraser
%   31      Tocantins
%   32      Vuoski & Neva
%
% Megan D Fowler, 2016-09-19
%   Updated 2016-10-04 to include additional rivers (29-32)
%   Created 2018-08-31 to get mask at 2 degree resolution 

coast = load('coast.mat');

nbasins = 32;   %Number of basins to consider

BasinNames={'YUKON'; 'MACKENZIE'; 'COLUMBIA'; 'NELSON'; 'MISSISSIPPI'; 'ST-LAWRENCE';'RHINE';'DANUBE';...
    'DNIEPR';'VOLGA';'OB';'YENISEI';'LENA';'AMUR';'ORINOCO';'PARANA';'AMAZONAS';'NIGER';'CONGO';'NILE';...
    'ZAMBEZI';'TIGRIS&EUPHRATES';'INDUS';'GANGES&BRAHMAPUTRA';'MEKONG';'HUANG-HE';'YANGTZE';'MURRAY&DARLING';...
    'DON'; 'FRASER'; 'TOCANTINS'; 'VUOKSI&NEVA'};

% % Read in CaMa data of lon and lat
% no_file = ['/Volumes/MyPassport/CaMa/output/noSP_cam3.5/outflwnoSP_cam3.5.nc'];
% lonDat  = ncread(no_file, 'lon');
% latDat  = ncread(no_file,'lat');
% camaRes = abs(latDat(2)-latDat(1)); %CaMa lat starts at the north pole

%Read in dimension sizes for files 
TRMMfile = '/Users/meganfowler/gp_fuse/obs/TRMM_1.9x2.5_199801-201312.nc';
lonDat = ncread(TRMMfile,'lon');
latDat = ncread(TRMMfile,'lat');
trmmRes = abs(latDat(2)-latDat(1)); 

%Convert longitude from 0:360 to -180:180
numlon = numel(lonDat);        %number of longitudes
delt_lon = lonDat(2)-lonDat(1);  %longitude resolution 
nlon = linspace(-180,(180-delt_lon),numlon)';   %Longitude array from -180 to just under 180

BasinMask=NaN(numel(lonDat),numel(latDat));

for ibasin=1:nbasins
    
    %Read in shape file for this particular basin
    filename=['/Volumes/MyPassport/CaMa/BasinDefinition-UN/Aqueduct_river_basins_',char(BasinNames(ibasin)),'/Aqueduct_river_basins_',char(BasinNames(ibasin))];
    M=m_shaperead(filename);       
    
    %Convert cell data of polygon to matlab array
    m      = cell2mat(M.ncst);
    bounds = M.mbr;
    %   Data is stored as an Nx2 array, with (:,1)=lon, and (:,2)=lat
    lonBasin = m(:,1);
    latBasin = m(:,2);
    
    %% Limit data to just one basin
    valcount=1;
    count=1;
    for i=1:numel(latBasin)
        
        %Check for missing data in the Basin latitudes
        if isnan(latBasin(i))==0
            [val index] = min(abs(latDat-latBasin(i)));  %CESM latitude closest to lat of basin shape file 
            ValidLat(valcount) = latDat(index);
            
            if i>=2 && ValidLat(valcount)~=ValidLat(valcount-1)
                NoDupeLat(count)=ValidLat(valcount-1); %Array of individual latitudes in the basin (in CaMa data)

                count=count+1;
            end

            valcount = valcount+1;
        end

    end

    for j=1:numel(NoDupeLat)
       k = find(latBasin>NoDupeLat(j)-trmmRes & latBasin<=NoDupeLat(j)+trmmRes); 

       lonsAtThisLat = lonBasin(k);
       lonMax        = max(lonsAtThisLat);
       lonMin        = min(lonsAtThisLat);

       [valMax indexMax] = min(abs(nlon-lonMax));
       [valMin indexMin] = min(abs(nlon-lonMin));

       LonMaxAtLat(j) = nlon(indexMax);
       LonMinAtLat(j) = nlon(indexMin);  

       klon = find(nlon>=LonMinAtLat(j) & nlon<=LonMaxAtLat(j));
       klat = find(latDat==NoDupeLat(j));

       BasinMask(klon,klat)= ibasin;
    end  
    
    clearvars ValidLat NoDupeLat LonMaxAtLat LonMinAtLat 
end

%Change longitudes from -180-180 to 0-360 
pnt0 = find(nlon==0);
basinMaskNew = [BasinMask(pnt0:end,:);BasinMask(1:pnt0-1,:)];

figure;
c=1:nbasins;
pcolor(lonDat,latDat,basinMaskNew')
%contourf(lonDat,latDat,basinMaskNew',c,'linecolor','none')
hold on;
plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k');
hold off;
axis([0 360 -90 90]);
title('TRMM resolution'); 
colorbar; 

% mekong  = nan(size(basinMaskNew));
% imekong = find(basinMaskNew==25);
% mekong(imekong)=1; 
% pcolor(lonDat,latDat,mekong');
% hold on;
% plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k');
% hold off;
% axis([0 360 -90 90]);


