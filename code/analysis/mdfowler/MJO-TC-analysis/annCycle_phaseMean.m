% Make plot of seasonal cycle of genesis 
%
% Run prep.m on individual phase's event dataset. Then run annual.m on
% that set, and save output: phs1_annCount = y;; clearvars -except
% phs1_annCount; 
% Repeat for each phase of the MJO. 
%
% Meg D. Fowler, 25 July 2019 

% Large array with genesis density in each phase
annCount_fullArr = nan(8,12);       % [phase, lat(y), lon(x)]              
annCount_fullArr(1,:) = phs1_annCount;
annCount_fullArr(2,:) = phs2_annCount;
annCount_fullArr(3,:) = phs3_annCount;
annCount_fullArr(4,:) = phs4_annCount;
annCount_fullArr(5,:) = phs5_annCount;
annCount_fullArr(6,:) = phs1_annCount;
annCount_fullArr(7,:) = phs7_annCount;
annCount_fullArr(8,:) = phs8_annCount;

%Take phase-mean
annCount_phsMean = squeeze(nanmean(annCount_fullArr,1));

%% Compare against obs and GPI
%  Observations can come from annualbesterr.m --> annCount_obs = y2;        *** annCount_obs ***
%  To get GPI average in each month, run GPI_PhaseMean_Bootstrapped.m 


%Limit GPI to same range and take mean/max --> *** GPI_monthly_phsMean ***
%   NOTE: GPI has a resolution of 2.5deg instead of 2, so not quite the
%   same lat/lon 

minLon = find(lon==(90));
maxLon = find(lon==(200));
minLat = find(lat==(30));
maxLat = find(lat==(0));
lat_GPI = lat(minLat:maxLat);
lon_GPI = lon(minLon:maxLon); 

GPI_region = GPI_monthly_phsMean(minLon:maxLon, minLat:maxLat,:);
GPI_monthly_mean = squeeze(nanmean(nanmean(GPI_region,1),2));


%% Plot differences
x = 1:12;   %Define months array

figure;
h1 = plot(x,annCount_obs, 'k','LineWidth',2,'Marker','s');
hold on;
h2 = plot(x, annCount_phsMean,'b','LineWidth',2,'Marker','s');
xlabel('Month');
ylabel('Normalized number of storms per month');
title('Annual cycle of cyclogenesis','FontSize',20)
xlim([0.5 12.5]);
grid on;
%Add GPI as new axis 
ax1 = gca; % current axes
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'YAxisLocation','right',...
    'Color','none');
ax2.YColor = 'red';
hold on;
%Plot GPI
h3 = plot(x,GPI_monthly_mean,'r','Marker','s','Parent',ax2,'linewidth',2);
ylabel('Monthly average GPI anomaly');
xlim([0.5 12.5]);
ylim([0 8])
legend([h1 h2 h3],'Best track', 'MIT model (Phase 1-8 average)', 'GPI anomaly');

%% Plot maps of seasonal amplitude and phase instead... 

% -- Define seasonal amplitude as max-min -- %
% Geographically, stored as longstore, latstore, and monthstore for each
%    event

% Modeled 

% -- Begin snippet from genpoints.m by Kerry Emanuel ---
params    %  Load parameters
%
clear x y vmask u
if exist('shearstore','var') == 0
    shape='circ';
    load('temp.mat')
end    
load('sorted')
[nn,m]=size(vstore);
%
if strcmp(bas,'GB')
    projection=gproject;
end    
%
nm=min(nn,nmax); % Determine number of points to plot
%
x=zeros(1,nm);
y=zeros(1,nm);
t=zeros(1,nm);
%
vmask=ceil(max(vmax-peakv,0)); % Creak mask of NaNs for storms that are insufficinetly intense
vmask=min(vmask,1);
vmask(vmask == 0)=NaN;
u=repmat(vmask,m,1);
%
xi=randperm(nn);   % Ranndomize squence of genesis points
longrand=u'.*longstore(xi,:);
latrand=u'.*latstore(xi,:);
timerand=u'.*monthstore(xi,:);   %MEG added this line
vrand=vnet(xi,:);
vmaxrand=vmax(xi);
%
[~,jmin]=min(max((startv-vrand),0),[],2); % Determine first points at which vnet exceeds startv
%
n=0;
for k=1:nm, 
   if vmaxrand(k) >= peakv  % Find first points at which vnet exceeds startv for storms with vmax > peakv
      n=n+1;
      x(n)=longrand(k,jmin(k));
      y(n)=latrand(k,jmin(k));
      t(n)=timerand(k,jmin(k));
   end   
end
% -- End part from K Emanuel --


% Sort into individual months... 
iJan   = find(t==1);
lonJan = x(iJan);
latJan = y(iJan); 
    
iFeb   = find(t==2);
lonFeb = x(iFeb);
latFeb = y(iFeb);

iMar   = find(t==3);
lonMar = x(iMar);
latMar = y(iMar); 

iApr   = find(t==4);
lonApr = x(iApr);
latApr = y(iApr); 

iMay   = find(t==5);
lonMay = x(iMay);
latMay = y(iMay);

iJun  = find(t==6);
lonJun = x(iJun);
latJun = y(iJun);

iJul = find(t==7); 
lonJul = x(iJul);
latJul = y(iJul); 

iAug = find(t==8);
lonAug = x(iAug);
latAug = y(iAug);

iSep = find(t==9); 
lonSep = x(iSep);
latSep = y(iSep); 

iOct = find(t==10);
lonOct = x(iOct);
latOct = y(iOct);

iNov = find(t==11);
lonNov = x(iNov);
latNov = y(iNov);

iDec = find(t==12);
lonDec = x(iDec);
latDec = y(iDec); 

% Plot  scatter with individual months in different colors 
coast = load('coast.mat');

figure; 
plot(coast.long, coast.lat, 'k','LineWidth',1.5);
oceancolor=[0.8 0.8 0.8];
set(gca,'color',oceancolor)
a=colormap(wmap);
for i=1:8
    a(i,:)=oceancolor;
end
colormap(a)
axis([90 200 0 35]);
hold on; 
%Plot each month
plot(lonJan,latJan,'.','Color',[1   0   0],'markersize',10); 
plot(lonFeb,latFeb,'.','Color',[1   0.5 0],'markersize',10); 
plot(lonMar,latMar,'.','Color',[1   1   0],'markersize',10); 
plot(lonApr,latApr,'.','Color',[0.5 1   0],'markersize',10); 
plot(lonMay,latMay,'.','Color',[0   1   0],'markersize',10); 
plot(lonJun,latJun,'.','Color',[0   1   0.5],'markersize',10); 
plot(lonJul,latJul,'.','Color',[0   1   1],'markersize',10); 
plot(lonAug,latAug,'.','Color',[0   0.5 1],'markersize',10); 
plot(lonSep,latSep,'.','Color',[0   0   1],'markersize',10); 
plot(lonOct,latOct,'.','Color',[0.5 0   1],'markersize',10); 
plot(lonNov,latNov,'.','Color',[1   0   1],'markersize',10); 
plot(lonDec,latDec,'.','Color',[1   0   0.5],'markersize',10); 
title('Genesis points by month','fontsize',20);
legend('Coast', 'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec');

%% Make gridded data and contour for May, Sep, and Dec 

% lonGrid = 100:2:180;
% latGrid = 0:2:30;
lonGrid = 100:5:180;
latGrid = 0:5:30;

for iLon=1:length(lonGrid)-1
    for iLat=1:length(latGrid)-1
       
        countMay(iLon,iLat) = length( find(lonMay>=lonGrid(iLon) & lonMay<lonGrid(iLon+1) & latMay>=latGrid(iLat) & latMay<latGrid(iLat+1)) );
        countSep(iLon,iLat) = length( find(lonSep>=lonGrid(iLon) & lonSep<lonGrid(iLon+1) & latSep>=latGrid(iLat) & latSep<latGrid(iLat+1)) );
        countDec(iLon,iLat) = length( find(lonDec>=lonGrid(iLon) & lonDec<lonGrid(iLon+1) & latDec>=latGrid(iLat) & latDec<latGrid(iLat+1)) );

    end
end

%Use midpoints for plotting
% lonPlot = 101:2:179;
% latPlot = 1:2:29;
lonPlot = 102.5:5:177.5;
latPlot = 2.5:5:27.5;

% try plotting 
coast = load('coast.mat');

figure; 
subplot(3,1,1)
contourf(lonPlot,latPlot,countMay',1:1:25,'linecolor','none');
hold on
plot(coast.long, coast.lat, 'k','LineWidth',1.5);
axis([100 180 0 30]);
cbar = colorbar();
caxis([1 25])
title('May Genesis in the MIT Model','fontsize',20);

%figure; 
subplot(3,1,2)
contourf(lonPlot,latPlot,countSep',1:1:45,'linecolor','none');
hold on
plot(coast.long, coast.lat, 'k','LineWidth',1.5);
axis([100 180 0 30]);
cbar = colorbar();
caxis([1 45])
title('September Genesis in the MIT Model','fontsize',20);

%figure; 
subplot(3,1,3)
contourf(lonPlot,latPlot,countDec',1:1:5,'linecolor','none');
hold on
plot(coast.long, coast.lat, 'k','LineWidth',1.5);
axis([100 180 0 30]);
cbar = colorbar();
caxis([1 5])
title('December Genesis in the MIT Model','fontsize',20);



% %Follow plotting approach of Kerry's
% % ---- MAY ---- %
% figure; 
% m_proj(projection,'long',[100 180],'lat',[0 30]);
% cmax=(max(max(countSep)));
% cint=cmax/15;
% [C1,h1]=m_contour(lonPlot,latPlot,countMay',0:1:cmax);
% %[C1,h1]=m_contour(x(kmin:kmax),y(imin:imax),genDen_phsMean(imin:imax,kmin:kmax),0:(0.1/15):0.1);   %Modified to have same colorbar
% colormap(wmap)
% if strcmp(wfill,'y')
%     hold on
%     h=m_pcolor(lonPlot,latPlot,countMay');
%     a=colormap(wmap);
%     for i=1:8,
%         a(i,:)=oceancolor;
%     end
%     colormap(a)
%     set(h,'EdgeAlpha',0,'FaceAlpha',wtrans);
%     set(h,'FaceColor','interp');
%     hold off
% end
% cb=colorbar;
% %set(cb,'Location','SouthOutside')    %Added by Meg
% set(cb,'fontweight',axisfontweight) 
% if strcmp(bas,'MT')
%      m_gshhs_l('patch',landcolor,'edgecolor','none');
% else     
%      m_coast('patch',landcolor,'edgecolor','none');
% end     
% if char(pstates) == char('y')
%     m_states(xmin, xmax, ymin, ymax,'color',stcolor)
% end 
% set(gca,'color',oceancolor)
% m_gridm('fontname',axisfont,'fontsize',20,'fontweight',axisfontweight,'linestyle',gridline)
% title('May Genesis in the MIT Model','fontsize',20);


%% Make similar plot for observations 

%Begin KE snippet from gendensitybesttrack.m 
params    %  Load parameters
%
clf('reset')
clear latscat longscat z2
%
pifac=acos(-1)/180;
coslati=1./(0.001+cos(pifac.*y));
%
bestproc   %  Process best tracks
%
[nn,m]=size(vbest);
%
[~,jmax]=min(vbest,[],2);
jmax=jmax-1;
%
[~,jmin]=min(max((startv-vbest),0),[],2);
%
yearbeg=min(bestyears);
yearend=max(bestyears);
%
n=0;
latscat=zeros(1,nn);
longscat=zeros(1,nn);
timescat = zeros(1,nn);  %ADDED by Meg 
for i=1:nn,
    if vmaxb(i) >= peakv 
      n=n+1;
      longscat(n)=longbest(i,jmin(i));
      latscat(n)=latbest(i,jmin(i));
      timescat(n) = monthbest(i,jmin(i)); % ADDED by Meg 
    end   
end
% -- End snippet 


% Sort into individual months... 
iJan_best   = find(timescat==1);
lonJan_best = longscat(iJan_best);
latJan_best = latscat(iJan_best); 
    
iFeb_best   = find(timescat==2);
lonFeb_best = longscat(iFeb_best);
latFeb_best = latscat(iFeb_best);

iMar_best   = find(timescat==3);
lonMar_best = longscat(iMar_best);
latMar_best = latscat(iMar_best); 

iApr_best   = find(timescat==4);
lonApr_best = longscat(iApr_best);
latApr_best = latscat(iApr_best); 

iMay_best   = find(timescat==5);
lonMay_best = longscat(iMay_best);
latMay_best = latscat(iMay_best);

iJun_best  = find(timescat==6);
lonJun_best = longscat(iJun_best);
latJun_best = latscat(iJun_best);

iJul_best = find(timescat==7); 
lonJul_best = longscat(iJul_best);
latJul_best = latscat(iJul_best); 

iAug_best = find(timescat==8);
lonAug_best = longscat(iAug_best);
latAug_best = latscat(iAug_best);

iSep_best = find(timescat==9); 
lonSep_best = longscat(iSep_best);
latSep_best = latscat(iSep_best); 

iOct_best = find(timescat==10);
lonOct_best = longscat(iOct_best);
latOct_best = latscat(iOct_best);

iNov_best = find(timescat==11);
lonNov_best = longscat(iNov_best);
latNov_best = latscat(iNov_best);

iDec_best = find(timescat==12);
lonDec_best = longscat(iDec_best);
latDec_best = latscat(iDec_best); 

figure; 
plot(coast.long, coast.lat, 'k','LineWidth',1.5);
oceancolor=[0.8 0.8 0.8];
set(gca,'color',oceancolor)
a=colormap(wmap);
for i=1:8
    a(i,:)=oceancolor;
end
colormap(a)
axis([90 200 0 35]);
hold on; 
%Plot each month
plot(lonJan_best,latJan_best,'.','Color',[1   0   0],'markersize',10); 
plot(lonFeb_best,latFeb_best,'.','Color',[1   0.5 0],'markersize',10); 
plot(lonMar_best,latMar_best,'.','Color',[1   1   0],'markersize',10); 
plot(lonApr_best,latApr_best,'.','Color',[0.5 1   0],'markersize',10); 
plot(lonMay_best,latMay_best,'.','Color',[0   1   0],'markersize',10); 
plot(lonJun_best,latJun_best,'.','Color',[0   1   0.5],'markersize',10); 
plot(lonJul_best,latJul_best,'.','Color',[0   1   1],'markersize',10); 
plot(lonAug_best,latAug_best,'.','Color',[0   0.5 1],'markersize',10); 
plot(lonSep_best,latSep_best,'.','Color',[0   0   1],'markersize',10); 
plot(lonOct_best,latOct_best,'.','Color',[0.5 0   1],'markersize',10); 
plot(lonNov_best,latNov_best,'.','Color',[1   0   1],'markersize',10); 
plot(lonDec_best,latDec_best,'.','Color',[1   0   0.5],'markersize',10); 
title('Genesis points by month [BestTrack]','fontsize',20);
legend('Coast', 'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec');


for iLon=1:length(lonGrid)-1
    for iLat=1:length(latGrid)-1
       
        countMay_best(iLon,iLat) = length( find(lonMay_best>=lonGrid(iLon) & lonMay_best<lonGrid(iLon+1) & latMay_best>=latGrid(iLat) & latMay_best<latGrid(iLat+1)) );
        countSep_best(iLon,iLat) = length( find(lonSep_best>=lonGrid(iLon) & lonSep_best<lonGrid(iLon+1) & latSep_best>=latGrid(iLat) & latSep_best<latGrid(iLat+1)) );
        countDec_best(iLon,iLat) = length( find(lonDec_best>=lonGrid(iLon) & lonDec_best<lonGrid(iLon+1) & latDec_best>=latGrid(iLat) & latDec_best<latGrid(iLat+1)) );

    end
end


% try plotting 
coast = load('coast.mat');

figure; 
subplot(3,1,1)
contourf(lonPlot,latPlot,countMay_best',1:1:5,'linecolor','none');
hold on
plot(coast.long, coast.lat, 'k','LineWidth',1.5);
axis([100 180 0 30]);
cbar = colorbar();
caxis([1 5])
title('May Genesis in Obs','fontsize',20);

%figure; 
subplot(3,1,2)
contourf(lonPlot,latPlot,countSep_best',1:1:15,'linecolor','none');
hold on
plot(coast.long, coast.lat, 'k','LineWidth',1.5);
axis([100 180 0 30]);
cbar = colorbar();
caxis([1 15])
title('September Genesis in Obs','fontsize',20);

%figure; 
subplot(3,1,3)
contourf(lonPlot,latPlot,countDec_best',1:1:5,'linecolor','none');
hold on
plot(coast.long, coast.lat, 'k','LineWidth',1.5);
axis([100 180 0 30]);
cbar = colorbar();
caxis([1 5])
title('December Genesis in obs','fontsize',20);

%% Make monthly maps of GPI as well... 

figure; 
subplot(3,1,1)
contourf(lon_GPI,lat_GPI,GPI_region(:,:,5)',1:1:10,'linecolor','none');
hold on
plot(coast.long, coast.lat, 'k','LineWidth',1.5);
axis([100 180 0 30]);
cbar = colorbar();
caxis([1 10])
title('May GPI','fontsize',20);
grid on; 

%figure; 
subplot(3,1,2)
contourf(lon_GPI,lat_GPI,GPI_region(:,:,9)',1:1:15,'linecolor','none');
hold on
plot(coast.long, coast.lat, 'k','LineWidth',1.5);
axis([100 180 0 30]);
cbar = colorbar();
caxis([1 15])
title('September Genesis in Obs','fontsize',20);
grid on;

%figure; 
subplot(3,1,3)
contourf(lon_GPI,lat_GPI,GPI_region(:,:,12)',1:1:5,'linecolor','none');
hold on
plot(coast.long, coast.lat, 'k','LineWidth',1.5);
axis([100 180 0 30]);
cbar = colorbar();
caxis([1 5])
title('December Genesis in obs','fontsize',20);
grid on;
