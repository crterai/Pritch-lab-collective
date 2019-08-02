% mapGPI_V2.m
%
% Create maps of GPI as calcualted for various basins in GPI_calc_v1.m.
% Plots are created for both no-SP and SP versions of CESM. This second
% version was created to use output of daily GPI computed via NCL remotely
% on Greenplanet. 
%
% Megan D. Fowler, 2017-03-24
%

%% Read in data
% CAMctrl_file   = '/Volumes/MyPassport/Data/TCs/GPI/CESM-CAM5-Control/CESM-CAM5-Control_AllTerms.nc';
% CAMexp_file    = '/Volumes/MyPassport/Data/TCs/GPI/CESM-CAM5-4xCO2/CESM-CAM5-4xCO2_AllTerms.nc';
% SPCAMctrl_file = '/Volumes/MyPassport/Data/TCs/GPI/CESM-SPCAM5-Control/CESM-SPCAM5-Control_AllTerms.nc';
% SPCAMexp_file  = '/Volumes/MyPassport/Data/TCs/GPI/CESM-SPCAM5-4xCO2/CESM-SPCAM5-4xCO2_AllTerms.nc';
 
CAMctrl_file   = '/Volumes/MyPassport/Data/TCs/GPI/CESM-CAM4-Control/CESM-CAM4-Control_AllTerms.nc';
CAMexp_file    = '/Volumes/MyPassport/Data/TCs/GPI/CESM-CAM4-4xCO2/CESM-CAM4-4xCO2_AllTerms.nc';
SPCAMctrl_file = '/Volumes/MyPassport/Data/TCs/GPI/CESM-SPCAM4-Control/CESM-SPCAM4-Control_AllTerms.nc';
SPCAMexp_file  = '/Volumes/MyPassport/Data/TCs/GPI/CESM-SPCAM4-4xCO2/CESM-SPCAM4-4xCO2_AllTerms.nc';

lon = ncread(CAMctrl_file,'lon');
lat = ncread(CAMctrl_file,'lat');

time_ctrl   = nctimenoleap(CAMctrl_file);
time_co2    = nctimenoleap(CAMexp_file);
time_SPctrl = nctimenoleap(SPCAMctrl_file);
time_SPco2  = nctimenoleap(SPCAMexp_file);

%Limit to exactly five years for CAM5 (Melissa runs)
% indices_timeCtrl   = find(time_ctrl >= datenum([1870 01 01 0 0 0])  & time_ctrl <= datenum([1875 01 01 0 0 0]));
% indices_timeCO2    = find(time_co2 >= datenum([2180 01 01 0 0 0])   & time_co2 <= datenum([2185 01 01 0 0 0]));
% indices_timeSPctrl = find(time_SPctrl >= datenum([1880 01 01 0 0 0])  & time_SPctrl <= datenum([1885 01 01 0 0 0]));
% indices_timeSPco2  = find(time_SPco2 >= datenum([2189 01 01 0 0 0])   & time_SPco2 <= datenum([2194 01 01 0 0 0]));

% %Limit to exactly five years for CAM4 (Mark runs)
indices_timeCtrl   = find(time_ctrl >= datenum([1880 01 01 0 0 0])  & time_ctrl <= datenum([1885 01 01 0 0 0]));
indices_timeCO2    = find(time_co2 >= datenum([2180 01 01 0 0 0])   & time_co2 <= datenum([2185 01 01 0 0 0]));
indices_timeSPctrl = find(time_SPctrl >= datenum([1880 01 01 0 0 0])  & time_SPctrl <= datenum([1885 01 01 0 0 0]));
indices_timeSPco2  = find(time_SPco2 >= datenum([2190 01 01 0 0 0])   & time_SPco2 <= datenum([2195 01 01 0 0 0]));

GPI_ctrl   = ncread(CAMctrl_file,'GPI');
GPI_co2    = ncread(CAMexp_file,'GPI');
GPI_SPctrl = ncread(SPCAMctrl_file,'GPI');
GPI_SPco2 = ncread(SPCAMexp_file,'GPI');

time_ctrl = time_ctrl(indices_timeCtrl);
time_co2  = time_co2(indices_timeCO2);
time_SPctrl = time_SPctrl(indices_timeSPctrl);
time_SPco2  = time_SPco2(indices_timeSPco2);

GPI_ctrl    = GPI_ctrl(:,:, indices_timeCtrl);
GPI_co2     = GPI_co2(:,:,indices_timeCO2);
GPI_SPctrl  = GPI_SPctrl(:,:,indices_timeSPctrl);
GPI_SPco2   = GPI_SPco2(:,:,indices_timeSPco2);

GPI_ctrl   = permute(GPI_ctrl, [3 2 1]);  %Time x lat x lon (just to match with below array set-up)
GPI_co2    = permute(GPI_co2, [3 2 1]);
GPI_SPctrl = permute(GPI_SPctrl, [3 2 1]);
GPI_SPco2  = permute(GPI_SPco2, [3 2 1]);

%% Compute monthly averages 
time_vecCtrl = datevec(double(time_ctrl));
yearCtrl     = time_vecCtrl(:,1);
monthCtrl    = time_vecCtrl(:,2);
JanCtrl = find(monthCtrl==2);
FebCtrl = find(monthCtrl==3);
MarCtrl = find(monthCtrl==4);
AprCtrl = find(monthCtrl==5);
MayCtrl = find(monthCtrl==6);
JunCtrl = find(monthCtrl==7);
JulCtrl = find(monthCtrl==8);
AugCtrl = find(monthCtrl==9);
SepCtrl = find(monthCtrl==10);
OctCtrl = find(monthCtrl==11);
NovCtrl = find(monthCtrl==12);
DecCtrl = find(monthCtrl==1);

GPIctrl_monthly  = [nanmean(GPI_ctrl(JanCtrl,:,:)); nanmean(GPI_ctrl(FebCtrl,:,:));...
    nanmean(GPI_ctrl(MarCtrl,:,:)); nanmean(GPI_ctrl(AprCtrl,:,:)); nanmean(GPI_ctrl(MayCtrl,:,:));nanmean(GPI_ctrl(JunCtrl,:,:));...
    nanmean(GPI_ctrl(JulCtrl,:,:)); nanmean(GPI_ctrl(AugCtrl,:,:)); nanmean(GPI_ctrl(SepCtrl,:,:));...
    nanmean(GPI_ctrl(OctCtrl,:,:)); nanmean(GPI_ctrl(NovCtrl,:,:)); nanmean(GPI_ctrl(DecCtrl,:,:))];


time_vecCO2  = datevec(double(time_co2));
yearCO2      = time_vecCO2(:,1);
monthCO2     = time_vecCO2(:,2);
JanCO2 = find(monthCO2==2);
FebCO2 = find(monthCO2==3);
MarCO2 = find(monthCO2==4);
AprCO2 = find(monthCO2==5);
MayCO2 = find(monthCO2==6);
JunCO2 = find(monthCO2==7);
JulCO2 = find(monthCO2==8);
AugCO2 = find(monthCO2==9);
SepCO2 = find(monthCO2==10);
OctCO2 = find(monthCO2==11);
NovCO2 = find(monthCO2==12);
DecCO2 = find(monthCO2==1);

GPIco2_monthly  = [nanmean(GPI_co2(JanCO2,:,:)); nanmean(GPI_co2(FebCO2,:,:));...
    nanmean(GPI_co2(MarCO2,:,:)); nanmean(GPI_co2(AprCO2,:,:)); nanmean(GPI_co2(MayCO2,:,:));nanmean(GPI_co2(JunCO2,:,:));...
    nanmean(GPI_co2(JulCO2,:,:)); nanmean(GPI_co2(AugCO2,:,:)); nanmean(GPI_co2(SepCO2,:,:));...
    nanmean(GPI_co2(OctCO2,:,:)); nanmean(GPI_co2(NovCO2,:,:)); nanmean(GPI_co2(DecCO2,:,:))];

time_vecSpCtrl  = datevec(double(time_SPctrl));
yearSpCtrl      = time_vecSpCtrl(:,1);
monthSPCtrl     = time_vecSpCtrl(:,2);
JanSpCtrl = find(monthSPCtrl==2);
FebSpCtrl = find(monthSPCtrl==3);
MarSpCtrl = find(monthSPCtrl==4);
AprSpCtrl = find(monthSPCtrl==5);
MaySpCtrl = find(monthSPCtrl==6);
JunSpCtrl = find(monthSPCtrl==7);
JulSpCtrl = find(monthSPCtrl==8);
AugSpCtrl = find(monthSPCtrl==9);
SepSpCtrl = find(monthSPCtrl==10);
OctSpCtrl = find(monthSPCtrl==11);
NovSpCtrl = find(monthSPCtrl==12);
DecSpCtrl = find(monthSPCtrl==1);

GPIspCtrl_monthly  = [nanmean(GPI_SPctrl(JanSpCtrl,:,:)); nanmean(GPI_SPctrl(FebSpCtrl,:,:));...
    nanmean(GPI_SPctrl(MarSpCtrl,:,:)); nanmean(GPI_SPctrl(AprSpCtrl,:,:)); nanmean(GPI_SPctrl(MaySpCtrl,:,:));nanmean(GPI_SPctrl(JunSpCtrl,:,:));...
    nanmean(GPI_SPctrl(JulSpCtrl,:,:)); nanmean(GPI_SPctrl(AugSpCtrl,:,:)); nanmean(GPI_SPctrl(SepSpCtrl,:,:));...
    nanmean(GPI_SPctrl(OctSpCtrl,:,:)); nanmean(GPI_SPctrl(NovSpCtrl,:,:)); nanmean(GPI_SPctrl(DecSpCtrl,:,:))];


time_vecSpCO2  = datevec(double(time_SPco2));
yearSpCO2      = time_vecSpCO2(:,1);
monthSPCO2     = time_vecSpCO2(:,2);
JanSpCO2 = find(monthSPCO2==2);
FebSpCO2 = find(monthSPCO2==3);
MarSpCO2 = find(monthSPCO2==4);
AprSpCO2 = find(monthSPCO2==5);
MaySpCO2 = find(monthSPCO2==6);
JunSpCO2 = find(monthSPCO2==7);
JulSpCO2 = find(monthSPCO2==8);
AugSpCO2 = find(monthSPCO2==9);
SepSpCO2 = find(monthSPCO2==10);
OctSpCO2 = find(monthSPCO2==11);
NovSpCO2 = find(monthSPCO2==12);
DecSpCO2 = find(monthSPCO2==1);

GPIspCO2_monthly  = [nanmean(GPI_SPco2(JanSpCO2,:,:)); nanmean(GPI_SPco2(FebSpCO2,:,:));...
    nanmean(GPI_SPco2(MarSpCO2,:,:)); nanmean(GPI_SPco2(AprSpCO2,:,:)); nanmean(GPI_SPco2(MaySpCO2,:,:));nanmean(GPI_SPco2(JunSpCO2,:,:));...
    nanmean(GPI_SPco2(JulSpCO2,:,:)); nanmean(GPI_SPco2(AugSpCO2,:,:)); nanmean(GPI_SPco2(SepSpCO2,:,:));...
    nanmean(GPI_SPco2(OctSpCO2,:,:)); nanmean(GPI_SPco2(NovSpCO2,:,:)); nanmean(GPI_SPco2(DecSpCO2,:,:))];

%% Compute annual maxima a different way 
%    This gives values that are WAY too large, so I really don't think that
%    this is the way to go on this one... 

% yrStartCtrl   = min(yearCtrl);
% yrStartCO2    = min(yearCO2);
% yrStartSPctrl = min(yearSpCtrl);
% yrStartSPco2  = min(yearSpCO2);
% 
% for iyr=1:5 
%     indicesCtrl_year   = find(yearCtrl==yrStartCtrl);
%     indicesCO2_year    = find(yearCO2==yrStartCO2);
%     indicesSPctrl_year = find(yearSpCtrl==yrStartSPctrl);
%     indicesSPco2_year  = find(yearSpCO2==yrStartSPco2);
%     
%     GPIctrl_year   = GPI_ctrl(indicesCtrl_year,:,:);
%     GPIco2_year    = GPI_co2(indicesCO2_year,:,:);
%     GPIspCtrl_year = GPI_SPctrl(indicesSPctrl_year,:,:);
%     GPIspCO2_year  = GPI_SPco2(indicesSPco2_year,:,:);
%     
%     AnnMaxGPI_ctrl(iyr,:,:)   = max(GPIctrl_year);
%     AnnMaxGPI_co2(iyr,:,:)    = max(GPIco2_year);
%     AnnMaxGPI_SPctrl(iyr,:,:) = max(GPIspCtrl_year);
%     AnnMaxGPI_SPco2(iyr,:,:)  = max(GPIspCO2_year);
%     
%     yrStartCtrl = yrStartCtrl+1;
%     yrStartCO2 = yrStartCO2+1;
%     yrStartSPctrl = yrStartSPctrl+1;
%     yrStartSPco2 = yrStartSPco2+1;
% end
% 
% avgAnnMax_ctrl   = squeeze(nanmean(AnnMaxGPI_ctrl));
% avgAnnMax_co2    = squeeze(nanmean(AnnMaxGPI_co2));
% avgAnnMax_SPctrl = squeeze(nanmean(AnnMaxGPI_SPctrl));
% avgAnnMax_SPco2  = squeeze(nanmean(AnnMaxGPI_SPco2));

%% Make Difference maps
coast = load('coast.mat');

% Use averages instead of annual max 
avgCtrl   = squeeze(nanmean(GPIctrl_monthly,1));
avgCO2    = squeeze(nanmean(GPIco2_monthly,1));
avgCtrlSP = squeeze(nanmean(GPIspCtrl_monthly,1));
avgCO2SP  = squeeze(nanmean(GPIspCO2_monthly,1));

%CAM5 v. CAM4 variables
% c5v4_control = avgCtrl_CAM5 - avgCtrl_CAM4; 
% c5v4_co2     = avgCO2_CAM5 - avgCO2_CAM4;
% c5v4_SPctrl  = avgCtrlSP_CAM5 - avgCtrlSP_CAM4;
% c5v4_SPco2   = avgCO2SP_CAM5 - avgCO2SP_CAM4;

%For a single version of CAM
timeChange_noSP = avgCO2-avgCtrl;
timeChange_SP   = avgCO2SP-avgCtrlSP;
SPchange_preind = avgCtrlSP-avgCtrl;
SPchange_future = avgCO2SP-avgCO2; 

figure;
subplot(2,1,1)
contourf(lon,lat,SPchange_preind,-7:0.07:7,'linecolor','none');
hold on;
plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k');
axis('xy','equal',[0 360 -45 45])
hold off;
%title('SPCAM5 - SPCAM4 [Control]');
%htitle=title('4xCO_2 - control [CESM1-CAM4]');
htitle=title('SP - noSP [control; CESM1-(SP)CAM4]');
colormap('jet');
h=colorbar('location','eastoutside');
caxis([-7 7])
set(h,'YTick',[-7 -5 -3 0 3 5 7]);

set(gca, 'FontName', 'Helvetica')
set(htitle,'FontName','Helvetica')
set(htitle,'FontSize',18)
set(gca, 'FontSize',16);
set(h, 'FontSize', 14);

subplot(2,1,2)
contourf(lon,lat,SPchange_future,-7:0.07:7,'linecolor','none');
hold on;
plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k');
axis('xy','equal',[0 360 -35 35])
hold off;
htitle=title('SP - noSP [4xCO_2; CESM1-(SP)CAM4]');
colormap('jet');
h2=colorbar('location','eastoutside');
caxis([-7 7])
set(h2,'YTick',[-7 -5 -3 0 3 5 7])

set(gca, 'FontName', 'Helvetica')
set(htitle,'FontName','Helvetica')
set(htitle,'FontSize',18)
set(gca, 'FontSize',16);
set(h2, 'FontSize', 14);

set(figure(1),'Position',[97 285 740 490])

%% Make maps
coast = load('coast.mat');

%%%%%% Conventional CESM %%%%%%%
% figure;
% subplot(2,1,1)
% contourf(lon,lat,squeeze(max(GPIctrl_monthly)),0:0.05:20,'linecolor','none');
% hold on;
% plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k');
% axis('xy','equal',[0 360 -45 45])
% hold off;
% title('Control GPI - CESM1-CAM4');
% colormap('jet');
% colorbar('location','eastoutside')
% 
% subplot(2,1,2)
% contourf(lon,lat,squeeze(max(GPIco2_monthly)),0:0.05:20,'linecolor','none');
% hold on;
% plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k');
% axis('xy','equal',[0 360 -35 35])
% hold off;
% title('4xCO_2 GPI - CESM1-CAM4');
% colormap('jet');
% colorbar('location','eastoutside')
% 
% %%%%%%% SP-CESM %%%%%%%
% figure;
% subplot(2,1,1)
% contourf(lon,lat,squeeze(max(GPIspCtrl_monthly)),0:0.05:20,'linecolor','none');
% hold on;
% plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k');
% axis('xy','equal',[0 360 -45 45])
% hold off;
% title('Control GPI - CESM1-SPCAM4');
% colormap('jet');
% colorbar('location','eastoutside')
% 
% subplot(2,1,2)
% contourf(lon,lat,squeeze(max(GPIspCO2_monthly)),0:0.05:20,'linecolor','none');
% hold on;
% plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k');
% axis('xy','equal',[0 360 -35 35])
% hold off;
% title('4xCO_2 GPI - CESM1-SPCAM4');
% colormap('jet');
% colorbar('location','eastoutside')

