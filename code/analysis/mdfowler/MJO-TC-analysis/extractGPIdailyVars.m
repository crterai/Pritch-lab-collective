% Adapted from extractDataForKerry to extract daily variables needed for
% computing GPI. 
%
% Meg D. Fowler, 2017-10-27
%

%% Set path to variable files and read in dimensions

load('/Users/meganfowler/Documents/MATLAB/TCs/SavedMatlabData/MJOindex&time_monthPhase');   % MJOindex_monthPhase, MJOtime_monthPhase

% ---------- Commented out below when computing climatology -------------
%load('/Users/meganfowler/Documents/MATLAB/TCs/SavedMatlabData/MJOindex&time_monthPhase-OMI');   % MJOindex_monthPhase, MJOtime_monthPhase
% -----------------------------------------------------------------------

SLPfiles = '/Volumes/MyPassport/Data/TCs/obs/DailySLP/slp_';

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

%for iPhase = 1:numel(retPhases)
for iPhase = 1:1
    %Pre-define arrays to save output in
    dailySST   = NaN(365,numel(lon),numel(lat));
    dailyT     = NaN(365,numel(lon),numel(lat),numel(levT));
    dailyQ     = NaN(365,numel(lon),numel(lat),numel(levT));
    dailyU     = NaN(365,numel(lon),numel(lat),numel(levU));
    dailyV     = NaN(365,numel(lon),numel(lat),numel(levU));
    dailySLP   = NaN(365,numel(lon),numel(lat));
    dates      = NaN(365);
    
    dayYr = 1.0;    %Day of the year
    
    for iMon = 1:numel(retMonths)
        %Use function to select which days to match with in obs record
        [getTimes]=selectLargeMJOevents(MJOindex_monthPhase,MJOtime_monthPhase,retMonths(iMon),retPhases(iPhase));
        
        ndays = numel(getTimes);    %Number of days in month
        for it = 1:numel(getTimes)
            retDate = find(ERAtime == getTimes(it));    %Date to retain
            
            % ---- For SLP
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

            % ---- For all other variables
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
        
        %Daily U and V values
        dailyU(dayYr:dayYr+(ndays-1),:,:,:) = Uday; %[ time lat lon lev]
        dailyV(dayYr:dayYr+(ndays-1),:,:,:) = Vday;
        %Re-arrange others so dimension order is same
        SSTday = permute(SSTday,[3 1 2]);
        Qday   = permute(Qday,[4 1 2 3]);
        Tday   = permute(Tday,[4 1 2 3]);
        SLPday = permute(SLPday,[3 1 2]);
        
        dailySST(dayYr:dayYr+(ndays-1),:,:) = SSTday; 
        dailyQ(dayYr:dayYr+(ndays-1),:,:,:) = Qday;
        dailyT(dayYr:dayYr+(ndays-1),:,:,:) = Tday; 
        dailySLP(dayYr:dayYr+(ndays-1),:,:) = SLPday;
        
        %Dates used in the averages and for U/V records
        dates(dayYr:dayYr+(ndays-1)) = getTimes;
        
        %Increment day of year to begin on the first of the next month
        dayYr = dayYr+ndays; 
        
        %Clearing variables necessary to avoid size errors on pre-existing arrays
        clearvars Uday Vday SSTday Qday Tday SLPday 
        
        %If verbose option is enabled, print status at end of each loop
        verbose=1;
        if verbose==1
            fprintf('Stored data for Phase %d month %d \n',retPhases(iPhase),iMon);
        end
    end
% Create netCDF files for phase 

       % Set missing values appropriately, rather than keeping as NaNs
       dailyU(isnan(dailyU))     = -32767;
       dailyV(isnan(dailyV))     = -32767;
       dailyT(isnan(dailyT))     = -32767;
       dailyQ(isnan(dailyQ))     = -32767;
       dailySST(isnan(dailySST)) = -32767;
       dailySLP(isnan(dailySLP)) = -32767;

       numPhase = num2str(retPhases(iPhase)); 
    
       filename = ['DailyGPIvars_MJOphase_',numPhase,'-RMM.nc'];

       nccreate(filename, 'lon',...
                'Dimensions',{'lon',numel(lon)});
            ncwrite(filename,'lon',lon);
            ncwriteatt(filename,'lon','long_name','degrees_east'); 

       nccreate(filename,'lat',...
                'Dimensions',{'lat',numel(lat)});
            ncwrite(filename,'lat',lat);
            ncwriteatt(filename,'lat','long_name','degrees_north');

       nccreate(filename,'levT',...
                'Dimensions',{'levT',numel(levT)});
            ncwrite(filename,'levT',levT);
            ncwriteatt(filename,'levT','long_name','Vertical pressure levels for T and Q');

       nccreate(filename,'levU',...
               'Dimensions',{'levU',numel(levU)});
         ncwrite(filename,'levU',levU);
         ncwriteatt(filename,'levU','long_name','Pressure levels for U and V');

        nccreate(filename,'T',...
                 'Dimensions',{'day',365,'lon',numel(lon),'lat',numel(lat),'levT',numel(levT)});
            ncwrite(filename,'T',dailyT);
            ncwriteatt(filename,'T','units','K');
            ncwriteatt(filename,'T','long_name','Temperature');
            ncwriteatt(filename,'T','FillValue', -32767);
            ncwriteatt(filename,'T','missing_value',-32767);

        nccreate(filename,'Q',...
                 'Dimensions',{'day',365,'lon',numel(lon),'lat',numel(lat),'levT',numel(levT)});
            ncwrite(filename,'Q',dailyQ);
            ncwriteatt(filename,'Q','units','kg/kg');
            ncwriteatt(filename,'Q','long_name','Specific humidity');
            ncwriteatt(filename,'Q','FillValue', -32767);
            ncwriteatt(filename,'Q','missing_value',-32767); 

        nccreate(filename,'SST',...
                 'Dimensions',{'day',365,'lon',numel(lon),'lat',numel(lat)});
             ncwrite(filename,'SST',dailySST);
             ncwriteatt(filename,'SST','units','K');
             ncwriteatt(filename,'SST','long_name','Sea surface temperature');
             ncwriteatt(filename,'SST','FillValue',-32767);
             ncwriteatt(filename,'SST','missing_value',-32767); 
             
        nccreate(filename,'SLP',...
                 'Dimensions',{'day',365,'lon',numel(lon),'lat',numel(lat)});
             ncwrite(filename,'SLP',dailySLP);
             ncwriteatt(filename,'SLP','units','Pa');
             ncwriteatt(filename,'SLP','long_name','Sea surface pressure');
             ncwriteatt(filename,'SLP','FillValue',-32767);
             ncwriteatt(filename,'SLP','missing_value',-32767); 

        nccreate(filename,'U',...
                 'Dimensions',{'day',365,'lon',numel(lon),'lat',numel(lat),'levU',numel(levU)});
             ncwrite(filename,'U',dailyU);
             ncwriteatt(filename,'U','units','m/s');
             ncwriteatt(filename,'U','long_name','U component of wind');
             ncwriteatt(filename,'U','FillValue',-32767);
             ncwriteatt(filename,'U','missing_value',-32767);

        nccreate(filename,'V',...
                 'Dimensions',{'day',365,'lon',numel(lon),'lat',numel(lat),'levU',numel(levU)});
             ncwrite(filename,'V',dailyV);
             ncwriteatt(filename,'V','units','m/s');
             ncwriteatt(filename,'V','long_name','U component of wind');
             ncwriteatt(filename,'V','FillValue',-32767);
             ncwriteatt(filename,'V','missing_value',-32767);     

         if verbose==1
             fprintf('File successfully created for phase of the MJO \n'); 
         end
    
end

%Save matlab variables in case I ever want to access them again, not in
%   netCDF format... 
%save('separatedObs_phases125','monthlySST','monthlyT','monthlyQ','dailyU','dailyV','dates');

%% Create climatology file 
% monthsArr = 1:12;
% daysArr   = [31 28 31 30 31 30 31 31 30 31 30 31];
% 
% vecTime = datevec(double(ERAtime));
% ERAmon = vecTime(:,2);
% ERAday = vecTime(:,3);
% 
% daycount  = 1;
% climSST   = NaN(365,numel(lon),numel(lat));
% climT     = NaN(365,numel(lon),numel(lat),numel(levT));
% climQ     = NaN(365,numel(lon),numel(lat),numel(levT));
% climU     = NaN(365,numel(lon),numel(lat),numel(levU));
% climV     = NaN(365,numel(lon),numel(lat),numel(levU));
% climSLP   = NaN(365,numel(lon),numel(lat));
% 
% for iMon = 1:numel(monthsArr)
%     tic
%     ndays = daysArr(iMon);    
%     for iday = 1:ndays 
%         iMatch = find(ERAmon==iMon & ERAday==iday); %Find where month/day line up
%         
%         for ii = 1:numel(iMatch)    %Loop over the 30 matches 
%             start = [1 1 1 iMatch(ii)];
%             stride = [1 1 1 1];
%             countU = [length(lon) length(lat) length(levU) length(iMatch(ii))]; 
%             countT = [length(lon) length(lat) length(levT) length(iMatch(ii))]; 
% 
%             vec = datevec(double(ERAtime(iMatch(ii))));
%             yr  = vec(1);
%             mon = vec(2);
%             day = vec(3);
%             
%             if mon<10 
%                slpFileName = [SLPfiles,'0',num2str(mon),'_',num2str(yr),'.nc'];
%             else
%                slpFileName = [SLPfiles,num2str(mon),'_',num2str(yr),'.nc']; 
%             end
%             
%             %Read in variables for selected date
%             SLP = ncread(slpFileName,'msl');    %[Pa]
%             SLPday(:,:,ii) = SLP(:,:,iday); 
%    
%             Uday(:,:,:,ii) = ncread(Ufile,'u',start,countU,stride);    %[m/s]
%             Vday(:,:,:,ii) = ncread(Vfile,'v',start,countU,stride);    %[m/s]
%             SSTday(:,:,ii) = ncread(SSTfile,'sst',[1 1 iMatch(ii)],[length(lon) length(lat) length(iMatch(ii))],[1 1 1]); %[K]
%             Qday(:,:,:,ii)   = ncread(Qfile,'q',start,countT,stride);  %[K]
%             Tday(:,:,:,ii)   = ncread(Tfile,'t',start,countT,stride);  %[kg/kg]
%             
%             %Replacing missing values with NaNs before averaging
%             Tday(Tday<=-32767)     = NaN;
%             Qday(Qday<=-32767)     = NaN;
%             SSTday(SSTday<=-32767) = NaN;
%             Uday(Uday<=-32767)     = NaN;
%             Vday(Vday<=-32767)     = NaN;
%             SLPday(SLPday<=-32767) = NaN;
% 
%         end
%         
%         % Compute climatology
%         climT(daycount,:,:,:) = squeeze(nanmean(Tday,4));
%         climQ(daycount,:,:,:) = squeeze(nanmean(Qday,4));
%         climU(daycount,:,:,:) = squeeze(nanmean(Uday,4));
%         climV(daycount,:,:,:) = squeeze(nanmean(Vday,4));
%         climSST(daycount,:,:) = squeeze(nanmean(SSTday,3));
%         climSLP(daycount,:,:) = squeeze(nanmean(SLPday,3)); 
%         
%         daycount = daycount+1;
%         
%         clearvars Tday Qday Uday Vday SSTday SLPday iMatch
%     end    
%     fprintf('Finished month %d \n',iMon);
%     toc
% end
% 
% % Set missing values appropriately, rather than keeping as NaNs
% climT(isnan(climT))     = -32767;
% climQ(isnan(climQ))     = -32767;
% climU(isnan(climU))     = -32767;
% climV(isnan(climV))     = -32767;
% climSST(isnan(climSST)) = -32767;
% climSLP(isnan(climSLP)) = -32767;
% 
% filename = ['DailyGPIvars_climatology_1983-2013.nc'];
% 
% nccreate(filename, 'lon',...
%         'Dimensions',{'lon',numel(lon)});
%     ncwrite(filename,'lon',lon);
%     ncwriteatt(filename,'lon','long_name','degrees_east'); 
% 
% nccreate(filename,'lat',...
%         'Dimensions',{'lat',numel(lat)});
%     ncwrite(filename,'lat',lat);
%     ncwriteatt(filename,'lat','long_name','degrees_north');
% 
% nccreate(filename,'levT',...
%         'Dimensions',{'levT',numel(levT)});
%     ncwrite(filename,'levT',levT);
%     ncwriteatt(filename,'levT','long_name','Vertical pressure levels for T and Q');
% 
% nccreate(filename,'levU',...
%        'Dimensions',{'levU',numel(levU)});
%  ncwrite(filename,'levU',levU);
%  ncwriteatt(filename,'levU','long_name','Pressure levels for U and V');
% 
% nccreate(filename,'T',...
%          'Dimensions',{'day',365,'lon',numel(lon),'lat',numel(lat),'levT',numel(levT)});
%     ncwrite(filename,'T',climT);
%     ncwriteatt(filename,'T','units','K');
%     ncwriteatt(filename,'T','long_name','Temperature');
%     ncwriteatt(filename,'T','FillValue', -32767);
%     ncwriteatt(filename,'T','missing_value',-32767);
% 
% nccreate(filename,'Q',...
%          'Dimensions',{'day',365,'lon',numel(lon),'lat',numel(lat),'levT',numel(levT)});
%     ncwrite(filename,'Q',climQ);
%     ncwriteatt(filename,'Q','units','kg/kg');
%     ncwriteatt(filename,'Q','long_name','Specific humidity');
%     ncwriteatt(filename,'Q','FillValue', -32767);
%     ncwriteatt(filename,'Q','missing_value',-32767); 
% 
% nccreate(filename,'SST',...
%          'Dimensions',{'day',365,'lon',numel(lon),'lat',numel(lat)});
%      ncwrite(filename,'SST',climSST);
%      ncwriteatt(filename,'SST','units','K');
%      ncwriteatt(filename,'SST','long_name','Sea surface temperature');
%      ncwriteatt(filename,'SST','FillValue',-32767);
%      ncwriteatt(filename,'SST','missing_value',-32767); 
% 
% nccreate(filename,'SLP',...
%          'Dimensions',{'day',365,'lon',numel(lon),'lat',numel(lat)});
%      ncwrite(filename,'SLP',climSLP);
%      ncwriteatt(filename,'SLP','units','Pa');
%      ncwriteatt(filename,'SLP','long_name','Sea surface pressure');
%      ncwriteatt(filename,'SLP','FillValue',-32767);
%      ncwriteatt(filename,'SLP','missing_value',-32767); 
% 
% nccreate(filename,'U',...
%          'Dimensions',{'day',365,'lon',numel(lon),'lat',numel(lat),'levU',numel(levU)});
%      ncwrite(filename,'U',climU);
%      ncwriteatt(filename,'U','units','m/s');
%      ncwriteatt(filename,'U','long_name','U component of wind');
%      ncwriteatt(filename,'U','FillValue',-32767);
%      ncwriteatt(filename,'U','missing_value',-32767);
% 
% nccreate(filename,'V',...
%          'Dimensions',{'day',365,'lon',numel(lon),'lat',numel(lat),'levU',numel(levU)});
%      ncwrite(filename,'V',climV);
%      ncwriteatt(filename,'V','units','m/s');
%      ncwriteatt(filename,'V','long_name','U component of wind');
%      ncwriteatt(filename,'V','FillValue',-32767);
%      ncwriteatt(filename,'V','missing_value',-32767);       
