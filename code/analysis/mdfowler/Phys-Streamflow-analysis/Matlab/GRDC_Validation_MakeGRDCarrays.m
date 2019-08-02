% This program takes the text files of data from GRDC, and extracts/saves
% the relevant information to a large matlab array that can be saved and
% accessed later. 
%
% Main point: The program creates a saved .mat array with summary
% information for each of the 163 basins supplied to us by GRDC. The format
% is as follows:
%   [stationID   River   lat   lon   catchmentArea   numberOfYearsInAvg   MeanAnnualMean   MeanAnnualMax   100yrDischarge (and now MeanAnnualMin)]
%
% Megan D Fowler, 2016-09-21

GRDCdir      = '/Volumes/MyPassport/Data/GRDC/data/'; %Directory where text files from GRDC are stored
basinIDarray = xlsread('~/Documents/Irvine/Flooding/GRDC-RequestedStations','A2:A164');  %read in basin IDs used in GRDC file names

nstations   = numel(basinIDarray);
GRDC_data    = cell(nstations,9);

%Constants for finding the 100-year discharge level
eulers = 0.57721;                                        %Euler's constant, used in computing Gumbel Distribution
K100   = (-sqrt(6)/pi)*(eulers + log(log(100/(100-1)))); %K value for 100-year discharge in Gumbel Distribution


for ibasin=1:nstations   
    ID = basinIDarray(ibasin);
    filename=[GRDCdir,num2str(ID),'.day'];
    fid = fopen(filename);
    
    line = fgetl(fid);
    %header_count=0;
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
        %header_count=header_count+1;    %Add one to counter tracking header lines
    end
        
    %Header lines not required since we've read them in now
    data      = textscan(fid, '%10s %*7c %8f %*1c %8f %*1c %4d','Delimiter',';');    
    date      = data{1};
    orig_data = data{2};        %m^3/s
    mod_data  = data{3};        %m^3/s
    flag      = data{4};       
    
    fclose(fid);                %Close open file
    %Assign initial basin-specific data (doesn't require calculations)
    GRDC_data(ibasin,1)={ID};
    GRDC_data(ibasin,2)={riverName};
    GRDC_data(ibasin,3)={latStation};
    GRDC_data(ibasin,4)={lonStation};
    GRDC_data(ibasin,5)={catchmentArea};
    
%     if ID==4214520
%         figure
%         days=1:numel(date);
%         plot(days,orig_data,'k');
%         xlabel('Time (days from 1)');
%         ylabel('Discharge (m^3/s)');
%         title(['Daily GRDC data from ', riverName]);
%     end
    
    %% Calculate data on station characteristics
    
    t = datetime(date,'InputFormat','yyyy-MM-dd');
    tvec = datevec(t);
    yr   = tvec(:,1);
    
    select = find(yr>=1970 & yr<=2005);
    
    yr       = yr(select);
    tvec     = tvec(select);    
    ryears   = yr(1):1:yr(end);     %Range of years being considered
    orig_sel = orig_data(select);
    mod_sel  = mod_data(select);
    flag_sel = flag(select);
    
    %Define arrays for this particular basin
    orig_annMean = NaN(numel(ryears),1);
    orig_annMax  = NaN(numel(ryears),1);
    orig_annMin  = NaN(numel(ryears),1);
    mod_annMean  = NaN(numel(ryears),1);
    mod_annMax   = NaN(numel(ryears),1);
    mod_annMin   = NaN(numel(ryears),1);
    
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
       
       orig_annMin(i)  = min(orig_year,[],'omitnan');
       mod_annMin(i)   = min(mod_year,[],'omitnan'); 
    end
    
    basinMean    = nanmean(orig_annMean);
    basinMeanMax = nanmean(orig_annMax);
    basinMeanMin = nanmean(orig_annMin);
        
    % Find 100-year discharge by using the Gumbel distribution
    qbar   = basinMeanMax;
    sigma  = std(orig_annMax,'omitnan');    
    q100   = qbar + (K100*sigma);
    
    %Assign rest of data to cell array
    GRDC_data(ibasin,6)={numel(ryears)};    %Number of years of data
    GRDC_data(ibasin,7)={basinMean};        %Mean discharge
    GRDC_data(ibasin,8)={basinMeanMax};     %Mean annual maximum discharge
    GRDC_data(ibasin,9)={q100};             %Discharge of 100-year flood
    GRDC_data(ibasin,10)={basinMeanMin};    %Mean annual minimum discharge
    
    % Create figures to save, showing annual mean/max discharge at each station
    figure;
    subplot(2,1,1);
    plot(ryears,mod_annMean);
    xlabel('Year');
    ylabel('Annual Mean Daily Discharge (m^3/s)');
    title(['Annual mean discharge for the ', riverName]);

    subplot(2,1,2);
    plot(ryears,mod_annMax);
    xlabel('Year');
    ylabel('Annual Mean Daily Discharge (m^3/s)');
    title(['Annual maximum discharge for the ', riverName]);
    
%    pdfname = ['~/Documents/MATLAB/FloodProject/Figures/GRDC-stationFigures/',riverName,num2str(ID)];
%    print(pdfname,'-dpng');
   
    close
end




