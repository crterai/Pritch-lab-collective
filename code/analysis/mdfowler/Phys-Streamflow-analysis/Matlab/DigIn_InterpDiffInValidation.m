% This program is to dig into why it makes such a difference to use linear
% interpolation vs. nearest neighbor when comparing GRDC observations with
% CaMa model output.
%
% Meg D Fowler, 2016-09-29
%

%% Read in data from both sources (model and obs)
%Read in modeled CaMa data
var     = 'outflw'; %Variable to analyze from annual max/mean .nc files
sp_file = ['/Volumes/MyPassport/CaMa/output/SP_cam3.5/',var,'SP_cam3.5.nc'];
no_file = ['/Volumes/MyPassport/CaMa/output/noSP_cam3.5/',var,'noSP_cam3.5.nc'];

lon  = ncread(sp_file, 'lon');
lat  = ncread(sp_file,'lat');
time = ncread(sp_file,'time');

sp_max  = ncread(sp_file,'ann_max'); %(m^3/s)
sp_mean = ncread(sp_file,'ann_mean');
no_max  = ncread(no_file,'ann_max');
no_mean = ncread(no_file,'ann_mean');

area = myarea(lon,lat);   %Modified grid cell area function from Gabe 

%Read in GRDC data
% s1ID = 4351900;    %Rio Grande
% s2ID = 2469260;    %Mekong
% s3ID = 6435060;    %Rhine

s1ID  = 4152050;     %Colorado
s2ID  = 2903430;     %Lena
s3ID  = 6742900;     %Danube

yrstart   = 1970;
yrend     = 1989;

[s1, s1_years, s1_lat, s1_lon, s1_ann_mean, s1_ann_max] = GRDC_get_yearlyData(s1ID,yrstart,yrend);
[s2, s2_years, s2_lat, s2_lon, s2_ann_mean, s2_ann_max] = GRDC_get_yearlyData(s2ID,yrstart,yrend);
[s3, s3_years, s3_lat, s3_lon, s3_ann_mean, s3_ann_max] = GRDC_get_yearlyData(s3ID,yrstart,yrend);

s1_lon = str2num(s1_lon);
s1_lat = str2num(s1_lat);
s2_lon = str2num(s2_lon);
s2_lat = str2num(s2_lat);
s3_lon = str2num(s3_lon);
s3_lat = str2num(s3_lat);

%% Find nearest neighbor, or interpolate to new values at the GRDC station 

%%%%%%%%%%     Nearest neighbor     %%%%%%%%%%%
[valLon, s1_iLon] = min(abs(lon - s1_lon));
[valLat, s1_iLat] = min(abs(lat - s1_lat));

[valLon, s2_iLon] = min(abs(lon - s2_lon));
[valLat, s2_iLat] = min(abs(lat - s2_lat));

[valLon, s3_iLon] = min(abs(lon - s3_lon));
[valLat, s3_iLat] = min(abs(lat - s3_lat));

s1_nearLon = lon(s1_iLon);            %Corresponding longitude and latitude in CaMa
s1_nearLat = lat(s1_iLat);

s2_nearLon = lon(s2_iLon);           
s2_nearLat = lat(s2_iLat);

s3_nearLon = lon(s3_iLon);            
s3_nearLat = lat(s3_iLat);

s1_nearMeanSP = squeeze(sp_mean(s1_iLon, s1_iLat, :));
s1_nearMeanNO = squeeze(no_mean(s1_iLon, s1_iLat, :));
s2_nearMeanSP = squeeze(sp_mean(s2_iLon, s2_iLat, :));
s2_nearMeanNO = squeeze(no_mean(s2_iLon, s2_iLat, :));
s3_nearMeanSP = squeeze(sp_mean(s3_iLon, s3_iLat, :));
s3_nearMeanNO = squeeze(no_mean(s3_iLon, s3_iLat, :));

s1_nearMaxSP = squeeze(sp_max(s1_iLon, s1_iLat, :));
s1_nearMaxNO = squeeze(no_max(s1_iLon, s1_iLat, :));
s2_nearMaxSP = squeeze(sp_max(s2_iLon, s2_iLat, :));
s2_nearMaxNO = squeeze(no_max(s2_iLon, s2_iLat, :));
s3_nearMaxSP = squeeze(sp_max(s3_iLon, s3_iLat, :));
s3_nearMaxNO = squeeze(no_max(s3_iLon, s3_iLat, :));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%   Linear interpolation   %%%%%%%%%%%
[X,Y] = meshgrid(lon,lat);
nyears = 20;

for t=1:nyears
    s1_SPmeanInterp(:,t) = interp2(X,Y,squeeze(sp_mean(:,:,t))',s1_lon,s1_lat,'linear');
    s1_NOmeanInterp(:,t) = interp2(X,Y,squeeze(no_mean(:,:,t))',s1_lon,s1_lat,'linear');
    s1_SPmaxInterp(:,t)  = interp2(X,Y,squeeze(sp_max(:,:,t))',s1_lon,s1_lat,'linear');
    s1_NOmaxInterp(:,t)  = interp2(X,Y,squeeze(no_max(:,:,t))',s1_lon,s1_lat,'linear');    

    s2_SPmeanInterp(:,t) = interp2(X,Y,squeeze(sp_mean(:,:,t))',s2_lon,s2_lat,'linear');
    s2_NOmeanInterp(:,t) = interp2(X,Y,squeeze(no_mean(:,:,t))',s2_lon,s2_lat,'linear');
    s2_SPmaxInterp(:,t)  = interp2(X,Y,squeeze(sp_max(:,:,t))',s2_lon,s2_lat,'linear');
    s2_NOmaxInterp(:,t)  = interp2(X,Y,squeeze(no_max(:,:,t))',s2_lon,s2_lat,'linear');    
       
    s3_SPmeanInterp(:,t) = interp2(X,Y,squeeze(sp_mean(:,:,t))',s3_lon,s3_lat,'linear');
    s3_NOmeanInterp(:,t) = interp2(X,Y,squeeze(no_mean(:,:,t))',s3_lon,s3_lat,'linear');
    s3_SPmaxInterp(:,t)  = interp2(X,Y,squeeze(sp_max(:,:,t))',s3_lon,s3_lat,'linear');
    s3_NOmaxInterp(:,t)  = interp2(X,Y,squeeze(no_max(:,:,t))',s3_lon,s3_lat,'linear');    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot all this stuff!
yr_range = 1:20;

figure;
plot(yr_range,s1_ann_mean,'k','LineWidth',1);
hold on;
plot(yr_range,s1_nearMeanSP,'b');
plot(yr_range,s1_SPmeanInterp,'--b');
plot(yr_range,s1_nearMeanNO,'r');
plot(yr_range,s1_NOmeanInterp,'--r');
hold off;
title(['Annual Mean Discharge - ', s1]);
legend('GRDC(1970-1989)','SP-nearest','SP-interpolated','noSP-nearest','noSP-interpolated');

figure;
plot(yr_range,s2_ann_mean,'k','LineWidth',1);
hold on;
plot(yr_range,s2_nearMeanSP,'b');
plot(yr_range,s2_SPmeanInterp,'--b');
plot(yr_range,s2_nearMeanNO,'r');
plot(yr_range,s2_NOmeanInterp,'--r');
hold off;
title(['Annual Mean Discharge - ', s2]);
legend('GRDC(1970-1989)','SP-nearest','SP-interpolated','noSP-nearest','noSP-interpolated');

figure;
plot(yr_range,s3_ann_mean,'k','LineWidth',1);
hold on;
plot(yr_range,s3_nearMeanSP,'b');
plot(yr_range,s3_SPmeanInterp,'--b');
plot(yr_range,s3_nearMeanNO,'r');
plot(yr_range,s3_NOmeanInterp,'--r');
hold off;
title(['Annual Mean Discharge - ', s3]);
legend('GRDC(1970-1989)','SP-nearest','SP-interpolated','noSP-nearest','noSP-interpolated');

%% Plot stations that DO have a difference and those that don't...
coast = load('coast.mat');

figure;
hold on;
contourf(lon,lat,squeeze(nanmean(no_mean,3))',100,'LineStyle','none');
plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k');
plot(-97.52, 25.9,'r*');    %Rio Grande
plot(6.11,51.84,'r*');      %Rhine
plot(28.7167, 45.2167, 'r*'); %Danube

plot(105.80,15.1167,'g*');  %Mekong
plot(-114.6319,32.7317,'g*');   %Colorado
plot(126.8, 72.37,'g*'); %Lena

hold off
axis('xy','equal',[-180 180 -90 90]);


