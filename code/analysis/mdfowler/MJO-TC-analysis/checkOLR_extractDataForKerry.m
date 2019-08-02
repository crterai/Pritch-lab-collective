% Check methodology of 'extractDataForKerry.m' by repeating the process
% using OLR data instead, and mapping out MJO-phase averages. 
%
% Meg D. Fowler, 2017-06-19

%% Set path to variable files and read in dimensions

load('/Volumes/MyPassport/Data/TCs/SavedMatlabData/MJOindex&time_monthPhase');   % MJOindex_monthPhase, MJOtime_monthPhase
load('/Volumes/MyPassport/Data/TCs/SavedMatlabData/OLRanom_new.mat');
load('/Volumes/MyPassport/Data/TCs/SavedMatlabData/timeForAnomalies_1983-2013_ERA-I.mat');


OLRfile    = '/Volumes/MyPassport/Data/TCs/obs/NOAA_interpolated_OLR_1983-2013.nc';
latNOAA = ncread(OLRfile, 'lat');
lonNOAA = ncread(OLRfile, 'lon');


%% Determine dates/phases to retreive 

retPhases = 1:8;    %Phases to retrieve
retMonths = 1:12;           %Months to retrieve

%Pre-define arrays to save output in
monthlyOLR = NaN(numel(retMonths),numel(retPhases),144,13);
dailyOLR   = NaN(365,numel(retPhases),144,13);
dates      = NaN(365,numel(retPhases));

for iPhase = 1:numel(retPhases)
    dayYr = 1.0;    %Day of the year
    for iMon = 1:numel(retMonths)
        %Use function to select which days to match with in obs record
        [getTimes]=selectLargeMJOevents(MJOindex_monthPhase,MJOtime_monthPhase,retMonths(iMon),retPhases(iPhase));
        
        ndays = numel(getTimes);    %Number of days in month
        for it = 1:numel(getTimes)
            retDate = find(time == getTimes(it));    %Date to retain
                        
            %Read in variables for selected date
            OLRmonth(it,:,:) = anomalyOLR(retDate,:,:);    %[m/s]
            OLRday(it,:,:)   = anomalyOLR(retDate,:,:);    %[m/s]             
        end
        
        %Average over the month and save result
        monthlyOLR(iMon,iPhase,:,:) = squeeze(nanmean(OLRmonth,1));
        
        %Daily U and V values
        dailyOLR(dayYr:dayYr+(ndays-1),iPhase,:,:) = OLRday; 
        
        %Dates used in the averages
        dates(dayYr:dayYr+(ndays-1),iPhase) = getTimes;
        
        %Increment day of year to begin on the first of the next month
        dayYr = dayYr+ndays; 
        
        %Clearing variables necessary to avoid size errors on pre-existing arrays
        clearvars OLRmonth OLRday
        
        %If verbose option is enabled, print status at end of each loop
        verbose=0;
        if verbose==1
            fprintf('Stored data for Phase %d month %d \n',retPhases(iPhase),iMon);
        end
    end
end

%Save matlab variables in case I ever want to access them again, not in
%   netCDF format... 
%save('separatedObs','monthlySST','monthlyT','monthlyQ','dailyU','dailyV','dates');

%% Check arrays... 
coast = load('coast.mat');
clearvars iMon 

% Check OLR progression per month

% c = -30:0.5:30;
% for iMon = 1:12
%     
%    figure; 
%    for iPhase = 1:8
%        OLRmonphs = squeeze(monthlyOLR(iMon,iPhase,:,:));
%        
%        subplot(8,1,iPhase);
%        contourf(lonNOAA,latNOAA,OLRmonphs',c,'LineColor','none');
%        title(['OLR anomaly for Month ', num2str(iMon),' & Phase ', num2str(iPhase)]);
%        hold on; 
%        plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k', 'LineWidth',2);
%        axis('xy','equal',[0 360 -20 20]);
%        %colormap('jet');
%        colorbar;
%        caxis([min(c) max(c)]);
% 
%    end
% end

wpacMonAvg = squeeze(nanmean(monthlyOLR(6:11,:,:,:),1));
c = -35:0.5:35;
figure;
for iPhase = 1:8
   subplot(8,1,iPhase);
   contourf(lonNOAA,latNOAA,squeeze(wpacMonAvg(iPhase,:,:))',c,'LineColor','none');
   title(['Jun-Nov OLR anomaly: Phase ', num2str(iPhase)]);
   hold on; 
   plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k', 'LineWidth',2);
   axis('xy','equal',[0 280 -20 20]);
   %colormap('jet');
   colorbar;
   caxis([min(c) max(c)]);

end

