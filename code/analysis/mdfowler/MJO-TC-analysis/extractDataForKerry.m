% Extract the needed variables for Kerry's downscaling model on days with
% an MJO index >=1. For pre-defined favorable/unfavorable MJO phases, some
% modifications have been made in order to ensure a more even distribution
% of the tails with respect to ENSO (Niño 3.4 index). 
%
% NEED TO HAVE GREENPLANET MOUNTED
%
% Meg D. Fowler, 2017-06-05
% Meg D. Fowler, 2017-08-25 - Modified to use OMI based MJO index and time
%   month-phase pairs. 
%

%% Set path to variable files and read in dimensions

%load('/Users/meganfowler/Documents/MATLAB/TCs/SavedMatlabData/MJOindex&time_monthPhase');   % MJOindex_monthPhase, MJOtime_monthPhase
load('/Users/meganfowler/Documents/MATLAB/TCs/SavedMatlabData/MJOindex&time_monthPhase-OMI');   % MJOindex_monthPhase, MJOtime_monthPhase

rootPath = '/Users/meganfowler/gp_fuse/TC/MJO_obs/ECMWF_ERA-Interim/2.5x2.5/'; %Path to Greenplanet ECMWF obs

Ufile    = [rootPath, 'dailyAvg/globalU_1983-2013_dailyAvg.nc'];
Vfile    = [rootPath, 'dailyAvg/globalV_1983-2013_dailyAvg.nc'];
SSTfile  = [rootPath, 'dailyAvg/globalSST_1983-2013_dailyAvg.nc'];
Tfile    = [rootPath, 'dailyAvg/globalT_1983-2013_dailyAvg.nc'];
Qfile    = [rootPath, 'dailyAvg/globalQ_1983-2013_dailyAvg.nc'];

%Read in time for all variables 
timeRaw       = ncread(Ufile,'time');        %Time is the same in U,V,Q,T, and SST files (confirmed)
startTime     = datenum([1900 01 01 0 0 0 ]);
ERAtime       = startTime + (timeRaw/24.0);      %Units of hours since 1900, but matlab operates on days

lon = ncread(SSTfile,'longitude');
lat = ncread(SSTfile,'latitude');
levU = ncread(Ufile, 'level');
levT = ncread(Tfile, 'level');

%% Determine dates/phases to retreive 

retPhases = 1:8;    %Phases to retrieve
retMonths = 1:12;           %Months to retrieve

%Pre-define arrays to save output in
monthlySST = NaN(numel(retMonths),numel(retPhases),numel(lon),numel(lat));
monthlyT   = NaN(numel(retMonths),numel(retPhases),numel(lon),numel(lat),numel(levT));
monthlyQ   = NaN(numel(retMonths),numel(retPhases), numel(lon),numel(lat),numel(levT));
dailyU     = NaN(365,numel(retPhases),numel(lon),numel(lat),numel(levU));
dailyV     = NaN(365,numel(retPhases),numel(lon),numel(lat),numel(levU));
dates      = NaN(365,numel(retPhases));

for iPhase = 1:numel(retPhases)
    dayYr = 1.0;    %Day of the year
    for iMon = 1:numel(retMonths)
        %Use function to select which days to match with in obs record
        [getTimes]=selectLargeMJOevents(MJOindex_monthPhase,MJOtime_monthPhase,retMonths(iMon),retPhases(iPhase));
        
        ndays = numel(getTimes);    %Number of days in month
        for it = 1:numel(getTimes)
            retDate = find(ERAtime == getTimes(it));    %Date to retain
            
            start = [1 1 1 retDate];
            stride = [1 1 1 1];
            countU = [length(lon) length(lat) length(levU) length(retDate)]; 
            countT = [length(lon) length(lat) length(levT) length(retDate)]; 
            
            %Read in variables for selected date
            Uday(it,:,:,:) = ncread(Ufile,'u',start,countU,stride);    %[m/s]
            Vday(it,:,:,:) = ncread(Vfile,'v',start,countU,stride);    %[m/s]
            SSTday(:,:,it) = ncread(SSTfile,'sst',[1 1 retDate],[length(lon) length(lat) length(retDate)],[1 1 1]); %[K]
            Qday(:,:,:,it)   = ncread(Qfile,'q',start,countT,stride);  %[K]
            Tday(:,:,:,it)   = ncread(Tfile,'t',start,countT,stride);  %[kg/kg]
            
            %Replacing missing values with NaNs before averaging
            Tday(Tday<=-32767)     = NaN;
            Qday(Qday<=-32767)     = NaN;
            SSTday(SSTday<=-32767) = NaN;
            Uday(Uday<=-32767)     = NaN;
            Vday(Vday<=-32767)     = NaN;
            
        end
        
        %Average over the month for SST, Q, and T; save result
        monthlySST(iMon,iPhase,:,:) = squeeze(nanmean(SSTday,3));
        monthlyQ(iMon,iPhase,:,:,:) = squeeze(nanmean(Qday,4));
        monthlyT(iMon,iPhase,:,:,:) = squeeze(nanmean(Tday,4));
        
        %Daily U and V values
        dailyU(dayYr:dayYr+(ndays-1),iPhase,:,:,:) = Uday; 
        dailyV(dayYr:dayYr+(ndays-1),iPhase,:,:,:) = Vday;
        
        %Dates used in the averages and for U/V records
        dates(dayYr:dayYr+(ndays-1),iPhase) = getTimes;
        
        %Increment day of year to begin on the first of the next month
        dayYr = dayYr+ndays; 
        
        %Clearing variables necessary to avoid size errors on pre-existing arrays
        clearvars Uday Vday SSTday Qday Tday
        
        %If verbose option is enabled, print status at end of each loop
        verbose=1;
        if verbose==1
            fprintf('Stored data for Phase %d month %d \n',retPhases(iPhase),iMon);
        end
    end
end

%Save matlab variables in case I ever want to access them again, not in
%   netCDF format... 
%save('separatedObs_phases125','monthlySST','monthlyT','monthlyQ','dailyU','dailyV','dates');

%% Write data out to netCDF file for each phase of the MJO 

% % Set missing values appropriately, rather than keeping as NaNs
% dailyU(isnan(dailyU))         = -32767;
% dailyV(isnan(dailyV))         = -32767;
% monthlyT(isnan(monthlyT))     = -32767;
% monthlyQ(isnan(monthlyQ))     = -32767;
% monthlySST(isnan(monthlySST)) = -32767;
% 
% verbose = 1;    %Print statement when done creating each netCDF file 
% 
% for iPhs=1:numel(retPhases)
%    numPhase = num2str(retPhases(iPhs)); 
%     
%    filename = ['OMI_annual_MJOphase',numPhase,'.nc'];
%    
%    nccreate(filename, 'lon',...
%             'Dimensions',{'lon',numel(lon)});
%         ncwrite(filename,'lon',lon);
%         ncwriteatt(filename,'lon','long_name','degrees_east'); 
%         
%    nccreate(filename,'lat',...
%             'Dimensions',{'lat',numel(lat)});
%         ncwrite(filename,'lat',lat);
%         ncwriteatt(filename,'lat','long_name','degrees_north');
% 
%    nccreate(filename,'levT',...
%             'Dimensions',{'levT',numel(levT)});
%         ncwrite(filename,'levT',levT);
%         ncwriteatt(filename,'levT','long_name','Vertical pressure levels for T and Q');
% 
%    nccreate(filename,'levU',...
%            'Dimensions',{'levU',numel(levU)});
%      ncwrite(filename,'levU',levU);
%      ncwriteatt(filename,'levU','long_name','Pressure levels for U and V');
%      
%     nccreate(filename,'T',...
%              'Dimensions',{'month',12,'lon',numel(lon),'lat',numel(lat),'levT',numel(levT)});
%         ncwrite(filename,'T',squeeze(monthlyT(:,iPhs,:,:,:)));
%         ncwriteatt(filename,'T','units','K');
%         ncwriteatt(filename,'T','long_name','Temperature');
%         ncwriteatt(filename,'T','FillValue', -32767);
%         ncwriteatt(filename,'T','missing_value',-32767);
%     
%     nccreate(filename,'Q',...
%              'Dimensions',{'month',12,'lon',numel(lon),'lat',numel(lat),'levT',numel(levT)});
%         ncwrite(filename,'Q',squeeze(monthlyQ(:,iPhs,:,:,:)));
%         ncwriteatt(filename,'Q','units','kg/kg');
%         ncwriteatt(filename,'Q','long_name','Specific humidity');
%         ncwriteatt(filename,'Q','FillValue', -32767);
%         ncwriteatt(filename,'Q','missing_value',-32767); 
%     
%     nccreate(filename,'SST',...
%              'Dimensions',{'month',12,'lon',numel(lon),'lat',numel(lat)});
%          ncwrite(filename,'SST',squeeze(monthlySST(:,iPhs,:,:)));
%          ncwriteatt(filename,'SST','units','K');
%          ncwriteatt(filename,'SST','long_name','Sea surface temperature');
%          ncwriteatt(filename,'SST','FillValue',-32767);
%          ncwriteatt(filename,'SST','missing_value',-32767); 
% 
%     nccreate(filename,'U',...
%              'Dimensions',{'day',365,'lon',numel(lon),'lat',numel(lat),'levU',numel(levU)});
%          ncwrite(filename,'U',squeeze(dailyU(:,iPhs,:,:,:)));
%          ncwriteatt(filename,'U','units','m/s');
%          ncwriteatt(filename,'U','long_name','U component of wind');
%          ncwriteatt(filename,'U','FillValue',-32767);
%          ncwriteatt(filename,'U','missing_value',-32767);
% 
%     nccreate(filename,'V',...
%              'Dimensions',{'day',365,'lon',numel(lon),'lat',numel(lat),'levU',numel(levU)});
%          ncwrite(filename,'V',squeeze(dailyV(:,iPhs,:,:,:)));
%          ncwriteatt(filename,'V','units','m/s');
%          ncwriteatt(filename,'V','long_name','U component of wind');
%          ncwriteatt(filename,'V','FillValue',-32767);
%          ncwriteatt(filename,'V','missing_value',-32767);     
%      
%      if verbose==1
%          fprintf('File successfully created for phase of the MJO \n'); 
%      end
%   
% end

%% Create file with surface pressure on same dates - used for GPI computation
SLPfiles = '/Volumes/MyPassport/Data/TCs/obs/DailySLP/slp_';
for iPhase = 1:numel(retPhases)
    dayYr = 1.0;    %Day of the year
    for iMon = 1:numel(retMonths)
        %Use function to select which days to match with in obs record
        [getTimes]=selectLargeMJOevents(MJOindex_monthPhase,MJOtime_monthPhase,retMonths(iMon),retPhases(iPhase));
        
        ndays = numel(getTimes);    %Number of days in month
        for it = 1:numel(getTimes)
            retDate = find(ERAtime == getTimes(it));    %Date to retain
            vec = datevec(double(ERAtime(retDate)));
            yr  = vec(1);
            mon = vec(2);
            day = vec(3);
            
            if mon<10 
               slpFileName = [SLPfiles,'0',num2str(mon),'_',num2str(yr),'.nc'];
            else
               slpFileName = [SLPfiles,num2str(mon),'_',num2str(yr),'.nc']; 
            end
            
            %Read in variables for selected date
            SLP = ncread(slpFileName,'msl');    %[Pa]
            SLPday(:,:,it) = SLP(:,:,day); 
        
            %Replacing missing values with NaNs before averaging
            SLPday(SLPday<=-32767)     = NaN;
            
        end
        
        %Average over the month for SST, Q, and T; save result
        monthlySST(iMon,iPhase,:,:) = squeeze(nanmean(SSTday,3));
        monthlyQ(iMon,iPhase,:,:,:) = squeeze(nanmean(Qday,4));
        monthlyT(iMon,iPhase,:,:,:) = squeeze(nanmean(Tday,4));
        
        %Daily U and V values
        dailyU(dayYr:dayYr+(ndays-1),iPhase,:,:,:) = Uday; 
        dailyV(dayYr:dayYr+(ndays-1),iPhase,:,:,:) = Vday;
        
        %Dates used in the averages and for U/V records
        dates(dayYr:dayYr+(ndays-1),iPhase) = getTimes;
        
        %Increment day of year to begin on the first of the next month
        dayYr = dayYr+ndays; 
        
        %Clearing variables necessary to avoid size errors on pre-existing arrays
        clearvars Uday Vday SSTday Qday Tday
        
        %If verbose option is enabled, print status at end of each loop
        verbose=1;
        if verbose==1
            fprintf('Stored data for Phase %d month %d \n',retPhases(iPhase),iMon);
        end
    end
end


%% Check data

monthlySST(monthlySST==-32767)=NaN;
monthlyT(monthlyT==-32767)=NaN;
monthlyQ(monthlyQ==-32767)=NaN;
dailyU(dailyU==-32767)=NaN;
dailyV(dailyV==-32767)=NaN;

cSST = 280:0.2:300;
cT   = 260:0.2:300;

coast = load('coast.mat');

for iMon=1:12
    figure;
    
    for iPhs=1:8
        SST = squeeze(monthlySST(iMon,iPhs,:,:));
        T   = squeeze(monthlyT(iMon,iPhs,:,:,end-1));   %975 mb
        Q   = squeeze(monthlyQ(iMon,iPhs,:,:,end-1));
        U   = squeeze(dailyU(iMon,iPhs,:,:,2)); %850 mb
        V   = squeeze(dailyV(iMon,iPhs,:,:,2));
        
    
       subplot(8,1,iPhs);
       contourf(lon,lat,T',cT,'LineColor','none');
       title(sprintf('T: Phase %d in month %d',iPhs,iMon));
%        hold on; 
%        plot(coast.long, coast.lat, 'k','LineWidth',2);
%        plot(coast.long+180, coast.lat, 'k','LineWidth',2);
       %axis('xy','equal',[0 180 -20 20]);
       colormap('jet');
       colorbar;
       %caxis([-15 15])
    end
end




