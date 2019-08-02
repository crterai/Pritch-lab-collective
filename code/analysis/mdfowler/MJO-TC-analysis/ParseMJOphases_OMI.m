% Compute MJO phase from OMI observations (1985-2013) and split into
% favorable/unfavorable phases for a given ocean basin. Based on previous
% script that used ERA-I data (ParseMJOphases_updated.m). 
% 
% Meg D. Fowler, 2017-08-24.
%


%% Determine MJO phase 
[mjoindex,mjophase,t] = OMI_obsPhaseIndex();

%Time from ERA-I RMMs 
timeVec = datevec(double(t));
yr      = timeVec(:,1);
mon     = timeVec(:,2);
day     = timeVec(:,3);

%Match the time period being analyzed to what we used for ERA-I RMMs 
matchStrt = find(t==datenum([1983 01 02]));
matchEnd  = find(t==datenum([2013 01 01])); 

t = t(matchStrt:matchEnd);
mjoindex = mjoindex(matchStrt:matchEnd);
mjophase = mjophase(matchStrt:matchEnd);

%Set threshold for determining MJO events 
threshold = 1;
%% Check that OLR looks reasonable in each phase
olrdebug=1;

if olrdebug==1
    %Use anomalous OLR

    %Read in OLR file to get lat/lon 
    OLRfile    = '/Volumes/MyPassport/Data/TCs/obs/NOAA_interpolated_OLR_1983-2013.nc';
    lat = ncread(OLRfile,'lat');
    lon = ncread(OLRfile,'lon');
    %Load coast file for plotting 
    coast = load('coast.mat');
    % Load anomalies of OLR and U for plotting (includes lat, which anomalies_ymean does not)
    load('OLRanom_new.mat');            %variable:"anomalyOLR"
    load('U200anomaly_withLat.mat');    %variable:"anomalyU200"
    load('U850anomaly_withLat.mat');    %variable:"anomalyU850"

    %Match OMI phase/amplitude with anomalies - remove first 120 days. 
    matchAnomPhase = mjophase(121:end);
    matchAnomIndex = mjoindex(121:end);
    matchAnomDays  = t(121:end);

    timeVec = datevec(matchAnomDays);
    mon     = timeVec(:,2);

    figure;
    c = -33:0.5:33;
    for iphase=1:8
       %it_mjo = find ((mon >= 11 | mon <= 5) & matchAnomIndex >= threshold & matchAnomPhase == iphase); %Winter
       it_mjo = find ((mon >= 6 | mon <= 11) & matchAnomIndex >= threshold & matchAnomPhase == iphase); % summer
       OLR_anom_map = squeeze(nanmean(anomalyOLR(it_mjo,:,:))); %Average of times matching condition
       subplot(8,1,iphase);
       contourf(lon,lat,OLR_anom_map',c,'LineColor','none');
       title(sprintf('June-Nov OLR anomaly: Phase %d',iphase));
%        if iphase==1
%           title(sprintf('anomalyOLR: Phase %d in ours [threshold = %d, winter]',iphase,threshold));
%        end
       hold on; 
       plot(coast.long, coast.lat, 'k','LineWidth',2);
       axis('xy','equal',[80 200 -20 20]); 
%       axis('xy','equal',[0 280 -20 20]); 
       %set(gca,'FontSize',12);
       colorbar;
       caxis([min(c) max(c)])

    end
end

%% Define MJO events

%Define MJO events as those above chosen threshold
iMJOevents    = find(mjoindex>=threshold); 
MJOeventIndex = mjoindex(iMJOevents);   %Extract indices of MJO events
MJOeventPhase = mjophase(iMJOevents);   %Extract phases of MJO events
MJOeventTime  = t(iMJOevents);       %Extract dates of MJO events

MJOeventTimeVec = datevec(double(MJOeventTime));
MJOeventMonth   = MJOeventTimeVec(:,2);

%% Separate months and MJO phases for basin's TC season

basin    = {'ATL',   'NE_Pac', 'NW_Pac',  'N_Ind',     'S_Ind',     'S_Pac'}; 
TCseason = {[7:11],  [6:10],   [6:11],    [4:6,10:12], [1:4,11,12], [1:4,12]}; %Cells with months of the TC season of each basin 
favHurr   = {[2],    [1,7,8],  [6,7],     [4,5],       [3,4],       [7]};      %Favorable phases for hurricanes per basin
unfavHurr = {[7],    [5],      [3,4],     [7,8],       [7,8],       [2,3,5]};  %Unfavorable phases for hurricanes per basin
    
% Create a cell array that's month x MJO phase, so that each month of the 
%   TC season has separate instances of each MJO phase. 
for iMon = 1:12  
    for iPhase=1:8  
        imonth_phase = find(MJOeventMonth==iMon & MJOeventPhase==iPhase); %Indices of where month and phase meet specification
        
        MJOindex_monthPhase{iMon,iPhase} = MJOeventIndex(imonth_phase);
        MJOtime_monthPhase{iMon,iPhase}  = MJOeventTime(imonth_phase); 
    end
end

%% Make plots of anomalous OLR propagation in each month

MJOdebug=0;

if MJOdebug==1
    %Use anomalous OLR

    %Read in OLR file to get lat/lon 
    OLRfile    = '/Volumes/MyPassport/Data/TCs/obs/NOAA_interpolated_OLR_1983-2013.nc';
    lat = ncread(OLRfile,'lat');
    lon = ncread(OLRfile,'lon');
    %Load coast file for plotting 
    coast = load('coast.mat');
    % Load anomalies of OLR and U for plotting (includes lat, which anomalies_ymean does not)
    load('OLRanom_new.mat');            %variable:"anomalyOLR"
    
    %OLR anomaly dates 
    matchAnomDays  = t(121:end);
    timeVec = datevec(matchAnomDays);
    mon     = timeVec(:,2);
    
    retPhases = 1:8;    %Phases to retrieve
    retMonths = 1:12;           %Months to retrieve
    
    for iMon=1:12
        for iPhase=1:8
            [getTimes]=selectLargeMJOevents(MJOindex_monthPhase,MJOtime_monthPhase,retMonths(iMon),retPhases(iPhase));

            ndays = numel(getTimes);    %Number of days in month
            for it = 1:ndays
                retDate = find(matchAnomDays == getTimes(it));    %Date to retain

                OLRphsPair(it,:,:) = anomalyOLR(retDate,:,:);            
            end

            monAvg(:,:,iMon,iPhase) = squeeze(nanmean(OLRphsPair,1));
        end
    end
    

    
    
    figure;
    for iphase=1:8
       subplot(8,1,iphase);
       contourf(lon,lat,squeeze(monAvg(:,:,12,iphase))','LineColor','none');
       title(sprintf('anomalyOLR: Phase %d in ours',iphase));
       if iphase==1
          title(sprintf('anomalyOLR: Phase %d in ours [threshold = %d, December]',iphase,threshold));
       end
       hold on; 
       plot(coast.long, coast.lat, 'k','LineWidth',2);
       axis('xy','equal',[0 180 -20 20]);
       colorbar;
       %caxis([-15 15])
    end
    
end


