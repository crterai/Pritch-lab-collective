% This program returns daily GRDC data for a specified station and year. 
%
% Inputs: station ID - corresponding to GRDC station ID 
%         year - specific year of interest, or if interested in all years, 'all'
%
% Meg D Fowler, 2016-10-06
% 

function [riverName, stationID, tvec, latStation, lonStation, orig_sel] = GRDC_get_dailyData(stationID,year)

%% Read in file
GRDCdir      = '/Volumes/MyPassport/Data/GRDC/data/';   %Directory where text files from GRDC are stored
filename     = [GRDCdir,num2str(stationID),'.day'];
fid = fopen(filename);

if fid==-1;    
    error('File not found. Check station ID. \n');
end

%Read in text in the header block
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

%Read in observed river discharges
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

yrRange = min(yr):max(yr);
for i=1:numel(yrRange);
   yrCount(i) = numel(find(yr==yrRange(i))); 
end

if year == 'all'
    orig_sel = orig_data; 
    
else    % If actual year is specified 
    
    select = find(yr==year);
    if numel(select)==0
        error('This station does not have data for the year specified.');
    end
    
    yr       = yr(select);
    tvec     = tvec(select);    
    orig_sel = orig_data(select);
    %mod_sel  = mod_data(select);

end


% yr       = yr(select);
% tvec     = tvec(select);    
% orig_sel = orig_data(select);
%mod_sel  = mod_data(select);

orig_sel(orig_sel==999)   = NaN;    %Replace missing values with NaN's
orig_sel(orig_sel==-999)  = NaN;    %Replace missing values with NaN's


end