% Compute MJO phase from ERA-I observations (1985-2013) and split into
% favorable/unfavorable phases for a given ocean basin. 
% 
% Meg D. Fowler, 2017-05-02 
%


%% Determine MJO phase 
load('anomalies_ymean_1983-2013_ERA-I.mat');        %Loads array "anomalies_ymean"
load('timeForAnomalies_1983-2013_ERA-I.mat');       %Loads time array

[RMM1,RMM2,RMM1components,RMM2components] = calculate_RMMs_from_anomalies_WHeofs (anomalies_ymean);
% RMM1var = std(RMM1);
% RMM2var = std(RMM2);
% RMM1    = RMM1/RMM1var;
% RMM2    = RMM2/RMM2var;

[mjoindex,mjophase] = RMMstophasespace(RMM1,RMM2);

iMJOevents = find(mjoindex>=1); 

MJOeventIndex = mjoindex(iMJOevents);   %Extract indices of MJO events
MJOeventPhase = mjophase(iMJOevents);   %Extract phases of MJO events
MJOeventTime  = time(iMJOevents);       %Extract dates of MJO events

%% Compare against WH observations... 
fid = fopen('rmmObs.74toRealtime.txt');
data = textscan(fid, '%d %d %d %f %f %d %f %s\n', 'Headerlines',2,'Delimiter' , '	','EndOfLine','\r\n');
fclose(fid);

yr_obs    = data{1};
mon_obs   = data{2};
day_obs   = data{3};
RMM1_obs  = data{4};
RMM2_obs  = data{5};
phase_obs = data{6};
amp_obs   = data{7};
timeObs   = datenum([double(yr_obs), double(mon_obs), double(day_obs)]);

% Use obs for now
% iMJOevents = find(amp_obs>=1); 
% 
% MJOeventIndex = amp_obs(iMJOevents);   %Extract indices of MJO events
% MJOeventPhase = phase_obs(iMJOevents); %Extract phases of MJO events
% MJOeventTime  = timeObs(iMJOevents);   %Extract dates of MJO events
% MJOeventMonth = mon_obs(iMJOevents);

% RMMs calculated in Mike's code...
timeVec = datevec(double(time));
yr      = timeVec(:,1);
mon     = timeVec(:,2);
day     = timeVec(:,3);

iMatch    = find(yr_obs==yr(1) & mon_obs==mon(1) & day_obs==day(1));
iMatchEnd = find(yr_obs==yr(end) & mon_obs==mon(end) & day_obs==day(end));
timeObs   = timeObs(iMatch:iMatchEnd);

RMM1_match = RMM1_obs(iMatch:iMatchEnd);
RMM2_match = RMM2_obs(iMatch:iMatchEnd); 
amp_match  = amp_obs(iMatch:iMatchEnd);
phs_match  = phase_obs(iMatch:iMatchEnd); 

OLRfile    = '/Volumes/MyPassport/Data/TCs/obs/NOAA_interpolated_OLR_1983-2013.nc';
lat = ncread(OLRfile,'lat');
lon = ncread(OLRfile,'lon');
coast = load('coast.mat');
threshold = 2;
load('OLRanom_new.mat');
load('U200anomaly_withLat.mat');
load('U850anomaly_withLat.mat');

event_debug = 1;
if (event_debug)    
    [x,y] = meshgrid (lon,timeObs(it));
    it = find (timeObs >= datenum([2007 8 1 0 0 0]) & timeObs <= datenum ([2008 5 28 0 0 0])); 
    figure;
    subplot (2,2,1); 
    plot(timeObs(it),amp_match(it),'k')
    hold on;
    plot(time(it),mjoindex(it),'r')
    %plot(time(it),1*ones(1,(3000-4)),'b'i,'LineWidth',2)
    datetick('x','mmm-dd-yyyy','keepticks');
    datetick ('x','mm/dd','keepticks'); 
    legend('WH Obs','ERA-I')
    grid on;
    subplot (2,2,2); 
    A = squeeze(nanmean(anomalyOLR(it,:,:),3)); % time x lon
    pcolor (x,y,A);
    shading flat; colorbar ('vert');  caxis ([-60 60]); 
    A = squeeze(nanmean(anomalyU200(it,:,:),3)); % time x lon
    hold on; 
    contour (x,y,A,[-8 -4 ],'k'); 
    contour (x,y,A,[4 8],'w'); 
    datetick ('y','mm/dd','keepticks'); 
    subplot (2,2,3); 
    %plot (RMM1(it),RMM2(it)); 
    plot (timeObs(it),mjophase(it));
    datetick ('x','mm/dd','keepticks'); xlabel ('Phase'); grid on;
    subplot (2,2,4);
    plot (timeObs(it),RMM1(it),'r'); hold on; plot (timeObs(it),RMM2(it),'b'); legend ({'RMM1','RMM2'}); 
    datetick ('x','mm/dd','keepticks'); 

    error ('debug pause'); 
end
% How does the propagation of the OLR appear for each MJO phase? 


%monClassic = find(mon==1 | mon==2 | mon==3 | mon==11 | mon==12);
%monClassic = find (mon >= 11 | mon <= 5); 
%anomalyOLR_winter = anomalyOLR(monClassic,:,:);
%anomalyU850_winter = anomalyU850(monClassic,:,:);
%anomalyU200_winter = anomalyU200(monClassic,:,:);
%mjoindex_winter   = mjoindex(monClassic);
%mjophase_winter   = mjophase(monClassic); 

%anomalyOLR_winter = anomalyOLR_winter./15.11623;    %Apply normalization factor?

%%
figure;
for iphase=1:8
%   it_mjo = find( mjoindex_winter>= threshold & mjophase_winter==iphase ); 
   it_mjo = find ((mon >= 11 | mon <= 5) & mjoindex' >= threshold & mjophase' == iphase); % Based on W-H RMMs.
   OLR_anom_map = squeeze(nanmean(anomalyOLR(it_mjo,:,:))); %Average of times matching condition
   subplot(8,1,iphase);
   contourf(lon,lat,OLR_anom_map','LineColor','none');
   title(sprintf('anomalyOLR: Phase %d in ours',iphase));
   if iphase==1
      title(sprintf('anomalyOLR: Phase %d in ours [threshold = %d, winter]',iphase,threshold));
   end
   hold on; 
   plot(coast.long, coast.lat, 'k','LineWidth',2);
   axis('xy','equal',[0 360 -20 20]);
   colorbar;
   %caxis([-15 15])

end


%% Set basin and favorable/unfavorable 

%Basin options: 
%   'ATL'   'NE_Pac'   'NW_Pac'   'N_Ind'   'S_Ind'   'S_Pac'
basin    = {'ATL', 'NE_Pac','NW_Pac',  'N_Ind',  'S_Ind',      'S_Pac'}; 
TCseason = {[7:11], [6:10], [6:11], [4:6,10:12],[1:4,11,12], [1:4,12]}; %Cells with months of the TC season of each basin 

favHurr   = {[2],[1,7,8],[6,7],[4,5],[3,4],[7,8]}; %favorable phases for hurricanes per basin
unfavHurr = {[7],[5],    [3,4],[7,8],[7,8],[2,3,5]}; %unfavorable phases for hurricanes per basin

for ibasin=1:numel(basin)       %Loop over basins
    basinName = basin{ibasin};  
    season    = TCseason{ibasin};
    
    % Create a cell array that's basin x month x MJO phase, so that for each
    %   basin, each month of the TC season has separate instances of each MJO
    %   phase. 
    for iMon = 1:numel(season)  %Loop over months in TC season
        for iPhase=1:8  %Loop over MJO phases
            month_phase{ibasin,iMon,iPhase} = find(MJOeventMonth==season(iMon) & MJOeventPhase==iPhase);        
        end
    end
    
    
    
end
