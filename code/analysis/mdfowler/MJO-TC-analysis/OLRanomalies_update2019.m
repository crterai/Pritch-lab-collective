% Look at OLR "anomalies" as OLR minus the phase 1-8 mean of OLR during
% those months. The previous version, "checkOLR_extractDataForKerry", used
% the traditional RMM approach which won't work necessarily for OMI and
% isn't consistent. That method computed the annual cycle with three
% harmonics maintained. That was then removed along with the mean of the
% prior 120 days. 
%
% Meg Fowler, 2019-05-13



load('/Volumes/MyPassport/Data/TCs/SavedMatlabData/MJOindex&time_monthPhase');   % MJOindex_monthPhase, MJOtime_monthPhase
load('/Volumes/MyPassport/Data/TCs/SavedMatlabData/timeForAnomalies_1983-2013_ERA-I.mat');

OLRfile    = '/Volumes/MyPassport/Data/TCs/obs/NOAA_interpolated_OLR_1983-2013.nc';
OLR = ncread(OLRfile,'olr');
latNOAA = ncread(OLRfile, 'lat');
lonNOAA = ncread(OLRfile, 'lon');

% To make sure times match up between OLR and MJO time series...
%    First day of MJO time series: 1983 - 05 - 02 - 0 - 0 - 0
%    OLR file time is in units: hours since 1800-01-01 00:00:0.0
timeOLR = ncread(OLRfile,'time');
baseTime = datenum(1800,1,1); 
% datevec(addtodate(baseTime,timeOLR(121),'hour')) --> When
%   timeOLR = 1983-05-02. So start there and go to end (end matches) 
timeOLR = timeOLR(121:end); 
OLR     = OLR(:,:,121:end);


%% Determine dates/phases to retreive 

retPhases = 1:8;    %Phases to retrieve
retMonths = 1:12;   %Months to retrieve

%Pre-define arrays to save output in
monthlyOLR = NaN(numel(retMonths),numel(retPhases),144,13);
dailyOLR   = NaN(numel(retPhases),144,13,365);
dates      = NaN(numel(retPhases),365);

for iPhase = 1:numel(retPhases)
    dayYr = 1.0;    %Day of the year
    for iMon = 1:numel(retMonths)
        %Use function to select which days to match with in obs record
        [getTimes]=selectLargeMJOevents(MJOindex_monthPhase,MJOtime_monthPhase,retMonths(iMon),retPhases(iPhase));
        
        ndays = numel(getTimes);    %Number of days in month
        for it = 1:numel(getTimes)
            retDate = find(time == getTimes(it));    %Date to retain
                        
            %Read in variables for selected date
            OLRmonth(:,:,it) = OLR(:,:,retDate);    %[m/s]
            OLRday(:,:,it)   = OLR(:,:,retDate);    %[m/s]             
        end
        
        %Average over the month and save result
        monthlyOLR(iMon,iPhase,:,:) = squeeze(nanmean(OLRmonth,3));
        
        %Daily U and V values
        dailyOLR(iPhase,:,:,dayYr:dayYr+(ndays-1)) = OLRday; 
        
        %Dates used in the averages
        dates(iPhase,dayYr:dayYr+(ndays-1)) = getTimes;
        
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

%% Average over all phases, Jun-Nov, and subtract to get 'anomaly' 

seasonalOLR = squeeze(nanmean(monthlyOLR(6:11,:,:,:)));
climOLR = squeeze(nanmean(seasonalOLR,1)); 

for iPhs = 1:8
   anomOLR(iPhs,:,:) = squeeze(seasonalOLR(iPhs,:,:))-climOLR;      
end

%% Make plots
coast = load('coast.mat');
clearvars iMon 

c = -40:0.5:40;
figure;
for iPhase = 1:8
   subplot(8,1,iPhase);
   contourf(lonNOAA,latNOAA,squeeze(anomOLR(iPhase,:,:))',c,'LineColor','none');
   title(['Jun-Nov OLR anomaly: Phase ', num2str(iPhase), ' (RMM)']);
   hold on; 
   plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k', 'LineWidth',2);
   axis('xy','equal',[0 280 -20 20]);
   colormap('jet');
   colorbar;
   caxis([min(c) max(c)]);

end

fig = gcf;
fig.PaperPositionMode = 'auto';
print('-painters','-depsc','~/Documents/Irvine/TCs/Figures/phaseAnomalies_OLR-RMM')


