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

coast = load('coast.mat');

nbasins = 32;   %Number of basins to consider

BasinNames={'YUKON'; 'MACKENZIE'; 'COLUMBIA'; 'NELSON'; 'MISSISSIPPI'; 'ST-LAWRENCE';'RHINE';'DANUBE';...
    'DNIEPR';'VOLGA';'OB';'YENISEI';'LENA';'AMUR';'ORINOCO';'PARANA';'AMAZONAS';'NIGER';'CONGO';'NILE';...
    'ZAMBEZI';'TIGRIS&EUPHRATES';'INDUS';'GANGES&BRAHMAPUTRA';'MEKONG';'HUANG-HE';'YANGTZE';'MURRAY&DARLING';...
    'DON'; 'FRASER'; 'TOCANTINS'; 'VUOKSI&NEVA'};


% Read in CaMa data of lon and lat
no_file = ['/Volumes/MyPassport/CaMa/output/noSP_cam3.5/outflwnoSP_cam3.5.nc'];
lonDat  = ncread(no_file, 'lon');
latDat  = ncread(no_file,'lat');
camaRes = abs(latDat(2)-latDat(1)); %CaMa lat starts at the north pole

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
            [val index] = min(abs(latDat-latBasin(i)));
            ValidLat(valcount) = latDat(index);
            
            if i>=2 & ValidLat(valcount)~=ValidLat(valcount-1)
                NoDupeLat(count)=ValidLat(valcount-1); %Array of individual latitudes in the basin (in CaMa data)

                count=count+1;
            end

            valcount = valcount+1;
        end

    end

    for j=1:numel(NoDupeLat)
       k = find(latBasin>NoDupeLat(j)-camaRes & latBasin<=NoDupeLat(j)+camaRes); 

       lonsAtThisLat = lonBasin(k);
       lonMax        = max(lonsAtThisLat);
       lonMin        = min(lonsAtThisLat);

       [valMax indexMax] = min(abs(lonDat-lonMax));
       [valMin indexMin] = min(abs(lonDat-lonMin));

       LonMaxAtLat(j) = lonDat(indexMax);
       LonMinAtLat(j) = lonDat(indexMin);  

       klon = find(lonDat>=LonMinAtLat(j) & lonDat<=LonMaxAtLat(j));
       klat = find(latDat==NoDupeLat(j));

       BasinMask(klon,klat)= ibasin;
    end  
    
    clearvars ValidLat NoDupeLat LonMaxAtLat LonMinAtLat 
end

figure;
c=1:nbasins;
contourf(lonDat,latDat,BasinMask',c)
hold on;
plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k');
hold off;
axis([-180 180 -90 90]);




