function [rivName, years, lat, lon, ann_mean, ann_max] = GRDC_get_yearlyData(stationID,yearstart,yearend)
% GRDC_get_yearlyData
% Given a GRDC station identifier, the function will (assuming it's a valid
% station) return a cell array of yearly annual max and mean flows. 
%
% INPUT: stationID - GRDC basin identifier used in file name
%        yearstart - First year the user wants retrieved (default = 1970)
%        yearend   - Last year the user wants retrieved  (default = 2005)
%
% Meg Fowler, 2016-09-29

%% Set defaults for year start & end 
if ~exist('yearstart','var')
    yearstart = 1970;
end
if ~exist('yearend','var')
    yearend = 2005;
end

%% Read in file
GRDCdir      = '/Volumes/MyPassport/Data/GRDC/data/';   %Directory where text files from GRDC are stored
filename     = [GRDCdir,num2str(stationID),'.day'];
fid = fopen(filename);

if fid==-1;    
    error('File not found. Check station ID. \n');
end

line = fgetl(fid);
while line(1)=='#'           %Marker for header lines in every text file (plus one additional line after that)
    if numel(line)>15        %Skip blank lines

        if line(3:7)=='River'   
            riverName = line(26:end);       %Save river name
            riverName(isspace(riverName)==1)='';    %Remove spaces from river name
            riverName(ismember(riverName,'/')) = '-';   %Replace slashes with dashes 
        elseif line(3:10)=='Latitude'
            latStation = line(26:end);      %Save station latitude, degrees
        elseif line(3:11)=='Longitude'      
            lonStation = line(26:end);      %Save station longitude, degrees
        elseif line(3:16)=='Catchment area'
            catchmentArea = line(26:end);   %Save station's catchment area, km^2 
        end

    end

    line = fgetl(fid);       %Read in next line
end

%Header lines not required since we've read them in now
data      = textscan(fid, '%10s %*7c %8f %*1c %8f %*1c %4d','Delimiter',';');    
date      = data{1};
orig_data = data{2};        %m^3/s
mod_data  = data{3};        %m^3/s
flag      = data{4};       

fclose(fid);                %Close open file
    
%% Calculate data on station characteristics
    
t = datetime(date,'InputFormat','yyyy-MM-dd');
tvec = datevec(t);
yr   = tvec(:,1);

select = find(yr>=yearstart & yr<=yearend);

yr       = yr(select);
tvec     = tvec(select);    
ryears   = yr(1):1:yr(end);     %Range of years being considered
orig_sel = orig_data(select);
mod_sel  = mod_data(select);
flag_sel = flag(select);
    
%Define arrays for this particular basin
orig_annMean = NaN(numel(ryears),1);
orig_annMax  = NaN(numel(ryears),1);
mod_annMean  = NaN(numel(ryears),1);
mod_annMax   = NaN(numel(ryears),1);
    
for i=1:numel(ryears)
   cyear = find(yr==ryears(i));      %Look at one year at a time - current year
   if yr(cyear)~=ryears(end)
       nyear = find(yr==ryears(i)+1);   %Indices of next year

       orig_year = orig_sel(cyear:(nyear-1));
       mod_year  = mod_sel(cyear:(nyear-1));
   elseif yr(cyear)==ryears(end)
       orig_year = orig_sel(cyear:end);
       mod_year  = mod_sel(cyear:end);           
   end

   orig_year(orig_year==-999) = NaN;    %Replace missing values with NaN's
   mod_year(mod_year==-999)   = NaN;
   orig_year(orig_year==999)  = NaN;    %Replace missing values with NaN's
   mod_year(mod_year==999)    = NaN;

   orig_annMean(i) = nanmean(orig_year);    %Mean annual discharge
   mod_annMean(i)  = nanmean(mod_year);     

   orig_annMax(i)  = max(orig_year,[],'omitnan');        %Annual maximum discharge
   mod_annMax(i)   = max(mod_year,[],'omitnan');     
end

rivName = riverName;
years = ryears;
lat   = latStation;
lon   = lonStation;
ann_mean = orig_annMean;
ann_max  = orig_annMax;

end


