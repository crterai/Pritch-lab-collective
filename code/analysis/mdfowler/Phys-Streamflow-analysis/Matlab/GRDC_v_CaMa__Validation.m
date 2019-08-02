% This program finds the annual mean, max, and 100-yr discharge at points
% corresponding (roughly) to those in the GRDC database. 
%
% Megan D Fowler, 2016-09-21
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
load('/Users/meganfowler/Documents/MATLAB/FloodProject/GRDC_data.mat');
stationID   = cell2mat(GRDC_data(:,1));
stationRiver = GRDC_data(:,2);
stationArea = cell2mat(GRDC_data(:,5));
stationlon  = cell2mat(GRDC_data(:,4));
stationlat  = cell2mat(GRDC_data(:,3));
stationMean = cell2mat(GRDC_data(:,7));
stationMax  = cell2mat(GRDC_data(:,8));
station100  = cell2mat(GRDC_data(:,9));

%% Match locations from model with those in GRDC
stationlon  = str2num(stationlon); %Originally stored as strings in GRDC text files
stationlat  = str2num(stationlat);
stationArea = str2num(stationArea);

%%%%%%Using interpolation to station locations%%%%%%
new_SPmean = permute(sp_mean, [3 2 1]);%Permute to be [time lat lon]
new_NOmean = permute(no_mean, [3 2 1]);
[X,Y] = meshgrid(lon,lat);

% 
% stationlon = sort(stationlon);
% stationlat = sort(stationlat);
% 
% 
% interpSP = ncl4matlab('linint2',lon,lat,new_SPmean,'True',stationlon,stationlat,0);

% interpSPmean = ncl4matlab('csa2l',lon,lat,new_SPmean,[4 4],stationlon,stationlat);

nyears = 20;

for t=1:nyears
    SPmeanInterp(:,t) = interp2(X,Y,squeeze(sp_mean(:,:,t))',stationlon,stationlat,'linear');
    NOmeanInterp(:,t) = interp2(X,Y,squeeze(no_mean(:,:,t))',stationlon,stationlat,'linear');
    SPmaxInterp(:,t)  = interp2(X,Y,squeeze(sp_max(:,:,t))',stationlon,stationlat,'linear');
    NOmaxInterp(:,t)  = interp2(X,Y,squeeze(no_max(:,:,t))',stationlon,stationlat,'linear');    

end
eulers = 0.57721;                                        %Euler's constant, used in computing Gumbel Distribution
K100   = -(sqrt(6)/pi)*(eulers + log(log(100/(100-1)))); %K value for 100-year discharge in Gumbel Distribution

spMax_pnt = squeeze(nanmean(SPmaxInterp,2));
noMax_pnt = squeeze(nanmean(NOmaxInterp,2));
spMean_pnt = squeeze(nanmean(SPmeanInterp,2));
noMean_pnt = squeeze(nanmean(NOmeanInterp,2));
qbarSP = spMean_pnt;
qbarNO = noMean_pnt;
sigmaSP = std(SPmaxInterp,0,2);
sigmaNO = std(NOmaxInterp,0,2);
q100_sp = qbarSP +(K100*sigmaSP);
q100_no = qbarNO + (K100*sigmaNO);



%% What about just the stations that are in H13? 

h13ID = [2903600; 3629000; 4125500; 4214270; 4152050; 4115200; 6742900;...
         5204251; 6978250; 2998400; 4122700; 2903430; 2469260; 4127800;...
         4122900; 4236010; 1234150; 2912600; 4123050; 2999910; 4122600;...
         4213681; 6435060; 4351900; 6970250; 4208400; 4116181; 6977100;...
         4120950; 2909150]; 
     
h13ID_exclusions = [2903600; 3629000; 4125500; 4214270; 6742900;...
         5204251; 6978250; 2998400; 4122700; 2903430; 2469260; 4127800;...
         4122900; 4236010; 1234150; 2912600; 4123050; 4122600;...
         4213681; 6435060; 4351900; 6970250; 4208400; 6977100;...
         4120950; 2909150]; %Removal of Columbia, Snake, Colorado, and Olenek. 
     
[boo, index] = ismember(h13ID,stationID);   %set this way on 2/1/18
%[boo, index] = ismember(h13ID_exclusions,stationID);

%What about in just the tropics? 
%index = find(stationlat>=-30 & stationlat<=30);

figure;
rno = corrcoef(stationMean(index),noMean_pnt(index));
rsp = corrcoef(stationMean(index),spMean_pnt(index));
loglog(stationMean(index),  noMean_pnt(index),'*r');
hold on;
loglog(stationMean(index), spMean_pnt(index), '*b');
hold off;
legend('noSP','SP');
title('Annual Mean Discharge');
ylabel('CaMa Discharge (m^3/s)');
xlabel('GRDC Discharge (m^3/s)');
ylim([1 10^6]);
xlim([1 10^6]);
hline=refline(1,0);
hline.Color='k';
text(1.5,10^4.7, ['noSP corrcoef: ', num2str(rno(1,2))]);
text(1.5,10^4.5, ['SP corrcoef:    ', num2str(rsp(1,2))]);

%Plot as separate
RMSE_noSP = sqrt(mean((noMean_pnt(index) - stationMean(index)).^2));  % Root Mean Squared Error
RMSE_SP   = sqrt(mean((spMean_pnt(index) - stationMean(index)).^2));  % Root Mean Squared Error

figure;
subplot(2,1,1)
rno = corrcoef(stationMean(index),noMean_pnt(index));
loglog(stationMean(index),  noMean_pnt(index),'*r');
hold on;
ylim([1 10^6]);
xlim([1 10^6]);
hline=refline(1,0);
hline.Color='k';
text(1.5,10^5.5, ['noSP corrcoef: ', num2str(rno(1,2))]);
text(1.5,10^5.0, ['noSP RMSE: ', num2str(RMSE_noSP)]);
title('noSP Validation')
subplot(2,1,2)
rsp = corrcoef(stationMean(index),spMean_pnt(index));
loglog(stationMean(index), spMean_pnt(index), '*b');
hold on;
ylim([1 10^6]);
xlim([1 10^6]);
hline=refline(1,0);
hline.Color='k';
text(1.5,10^5.5, ['SP corrcoef:    ', num2str(rsp(1,2))]);
text(1.5,10^5.0, ['SP RMSE: ', num2str(RMSE_SP)]);

title('SP Validation')






% 
% figure;
% rno = corrcoef(stationMax(index),noMax_pnt(index));
% rsp = corrcoef(stationMax(index),spMax_pnt(index));
% loglog(stationMax(index),  noMax_pnt(index),'*r');
% hold on;
% loglog(stationMax(index), spMax_pnt(index), '*b');
% hold off;
% legend('noSP','SP');
% title('Annual Maximum Discharge');
% ylabel('CaMa Discharge (m^3/s)');
% xlabel('GRDC Discharge (m^3/s)');
% ylim([1 10^6]);
% xlim([1 10^6]);
% hline=refline(1,0);
% hline.Color='k';
% text(1.5,10^5.7, ['noSP corrcoef: ', num2str(rno(1,2))]);
% text(1.5,10^5.5, ['SP corrcoef:    ', num2str(rsp(1,2))]);
% 
% 
% figure;
% rno = corrcoef(station100(index),q100_no(index));
% rsp = corrcoef(station100(index),q100_sp(index));
% loglog(station100(index),  q100_no(index),'*r');
% hold on;
% loglog(station100(index), q100_sp(index), '*b');
% hold off;
% legend('noSP','SP');
% title('100-Year Discharge');
% ylabel('CaMa Discharge (m^3/s)');
% xlabel('GRDC Discharge (m^3/s)');
% ylim([1 10^6]);
% xlim([1 10^6]);
% hline=refline(1,0);
% hline.Color='k';
% text(1.5,10^5.7, ['noSP corrcoef: ', num2str(rno(1,2))]);
% text(1.5,10^5.5, ['SP corrcoef:    ', num2str(rsp(1,2))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% What about rivers that are heavily vs. weakly managed? 

% % heavyID = [2909150;2909152;4127501;4127800;4115200;4115201;6977100];
% % weakID  = [2569005;2969010;2469050;2469072;2969090;2969095;2969100;2969101;...
% %             2469260;2969430;3629000;3629001;3206720];
%         
% heavyID = [1291100;1891500;1234150;1734500;6978250;4127501;4127800;...
%            4115200;4115201;4213720;3265601;3649900;3649950;2909150;2909152;6972430];
% weakID = [1147010;2469260;6335020;6435060;3629000;3629001;3206720;4103200;4103550];
% 
% 
% [boolh, iHeavy] = ismember(heavyID,stationID);
% [boolw, iWeak] = ismember(weakID,stationID);
% 
% figure;
% rno = corrcoef(stationMax(iHeavy),noMax_pnt(iHeavy));
% rsp = corrcoef(stationMax(iHeavy),spMax_pnt(iHeavy));
% loglog(stationMax(iHeavy),  noMax_pnt(iHeavy),'*r');
% hold on;
% loglog(stationMax(iHeavy), spMax_pnt(iHeavy), '*b');
% hold off;
% legend('noSP','SP');
% title('Annual Max Discharge - Heavily Managed');
% ylabel('CaMa Discharge (m^3/s)');
% xlabel('GRDC Discharge (m^3/s)');
% % ylim([1 10^5]);
% % xlim([1 10^5]);
% hline=refline(1,0);
% hline.Color='k';
% text(1.5,10^4.7, ['noSP corrcoef: ', num2str(rno(1,2))]);
% text(1.5,10^4.5, ['SP corrcoef:    ', num2str(rsp(1,2))]);
% 
% figure;
% rno = corrcoef(stationMax(iWeak),noMax_pnt(iWeak));
% rsp = corrcoef(stationMax(iWeak),spMax_pnt(iWeak));
% loglog(stationMax(iWeak),  noMax_pnt(iWeak),'*r');
% hold on;
% loglog(stationMax(iWeak), spMax_pnt(iWeak), '*b');
% hold off;
% legend('noSP','SP');
% title('Annual Max Discharge - Lightly Managed');
% ylabel('CaMa Discharge (m^3/s)');
% xlabel('GRDC Discharge (m^3/s)');
% % ylim([1 10^5]);
% % xlim([1 10^5]);
% hline=refline(1,0);
% hline.Color='k';
% text(1.5,10^4.7, ['noSP corrcoef: ', num2str(rno(1,2))]);
% text(1.5,10^4.5, ['SP corrcoef:    ', num2str(rsp(1,2))]);

% nstations  = numel(stationID);       %Number of stations in record
% 
% radEarth = 6371.0;  %radius of Earth in km
% 
% %Loop over stations and find closest CaMa grid cells 
% for s=1:nstations
%     [valLon indexLon] = min(abs(lon - stationlon(s)));
%     [valLat indexLat] = min(abs(lat - stationlat(s)));
% 
%     cLon(s) = lon(indexLon);            %Corresponding longitude and latitude in CaMa
%     cLat(s) = lat(indexLat);
%     
%     cellArea(s) = area(indexLon,indexLat);
%     
%     dlon = abs(stationlon(s)-cLon(s));  %Convert from degrees to radians
%     dlat = abs(stationlat(s)-cLat(s));
%     
%     %Compute great circle distance between the two points using the
%     %haversine formula 
%     angleBtw = 2*asin(sqrt( (sin(dlat/2)).^2 + cos(cLat(s))*cos(stationlat(s))* (sin(dlon/2)).^2 ));
%     grtcircledist(s) = (angleBtw*(pi/180))*radEarth;
%     
% end
% % figure;
% % coast = load('coast.mat');
% % hold on;
% % plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k');
% % plot(stationlon,stationlat,'.r')
% % plot(cLon,cLat,'.b')
% % hold off
% % axis('xy','equal',[-180 180 -90 90]);
% fprintf('The mean distance between the two data points is %f km \n', nanmean(grtcircledist));
% 
% %% Compute annual mean, mean annual max, and 100-year discharge for each point
% 
% eulers = 0.57721;                                        %Euler's constant, used in computing Gumbel Distribution
% K100   = -(sqrt(6)/pi)*(eulers + log(log(100/(100-1)))); %K value for 100-year discharge in Gumbel Distribution
% 
% SP_data    = cell(nstations,9);
% noSP_data  = cell(nstations,9);
% 
% ryears = 1:20;
% 
% for s=1:nstations
%     ilon = find(lon==cLon(s));
%     ilat = find(lat==cLat(s));
%     
%     spMax_pnt  = squeeze(sp_max(ilon,ilat,:));     %[m^3/s]
%     noMax_pnt  = squeeze(no_max(ilon,ilat,:));
%     spMean_pnt = squeeze(sp_mean(ilon,ilat,:));
%     noMean_pnt = squeeze(no_mean(ilon,ilat,:));
%     
%     pntMeanSP(s) = squeeze(nanmean(spMean_pnt));
%     pntMeanNo(s) = squeeze(nanmean(noMean_pnt));
%     meanMaxSP(s) = squeeze(nanmean(spMax_pnt));
%     meanMaxNo(s) = squeeze(nanmean(noMax_pnt));
%     
%     %Compute 100-year discharge
%     qbarSP = meanMaxSP(s);
%     qbarNo = meanMaxNo(s);
%     sigmaSP = std(spMax_pnt);
%     sigmaNo = std(noMax_pnt);
%     q100_sp(s) = qbarSP + (K100*sigmaSP);
%     q100_no(s) = qbarNo + (K100*sigmaNo);    
%     
%     %% save SP/noSP to cell arrays, similarly to for GRDC data
% 
%     SP_data(s,1) = {stationID(s)};
%     SP_data(s,2) = {stationRiver(s)};
%     SP_data(s,3) = {cLat(s)};
%     SP_data(s,4) = {cLon(s)};
%     SP_data(s,5) = {cellArea(s)};
%     SP_data(s,6) = {20};
%     SP_data(s,7) = {pntMeanSP(s)};
%     SP_data(s,8) = {meanMaxSP(s)};
%     SP_data(s,9) = {q100_sp(s)};
% 
%     noSP_data(s,1) = {stationID(s)};
%     noSP_data(s,2) = {stationRiver(s)};
%     noSP_data(s,3) = {cLat(s)};
%     noSP_data(s,4) = {cLon(s)};
%     noSP_data(s,5) = {cellArea(s)};
%     noSP_data(s,6) = {20};
%     noSP_data(s,7) = {pntMeanNo(s)};
%     noSP_data(s,8) = {meanMaxNo(s)};
%     noSP_data(s,9) = {q100_no(s)};
%     
%     % Create figures to save, showing annual mean/max discharge at each station
% %     figure;
% %     subplot(2,1,1);
% %     plot(ryears,noMean_pnt,'r');
% %     hold on;
% %     plot(ryears,spMean_pnt,'b');
% %     hold off;
% %     legend('noSP','SP');
% %     xlabel('Year');
% %     ylabel('Annual Mean Daily Discharge (m^3/s)');
% %     title(['Annual mean discharge for the ', stationRiver{s}]);
% % 
% %     subplot(2,1,2);
% %     plot(ryears,noMax_pnt,'r');
% %     hold on;
% %     plot(ryears,spMax_pnt,'b');
% %     hold off
% %     legend('noSP','SP')
% %     xlabel('Year');
% %     ylabel('Annual Mean Daily Discharge (m^3/s)');
% %     title(['Annual maximum discharge for the ', stationRiver{s}]);
% %     
% %     pdfname = ['~/Documents/MATLAB/FloodProject/Figures/CaMa-StationFigures/',stationRiver{s},num2str(stationID(s))];
% %     print(pdfname,'-dpng');
% %    
% %     close
% end
% 
% % Normalize by area of catchment/grid cell (converts units to m/s instead
% % of m^3/s) 
% 
% %%%%% This seems to make agreement between the two way worse. %%%%%%%%%%%
% 
% stationArea = stationArea*(10^6); %Convert from km^2 to m^2
% cellArea    = cellArea*(10^6);    %Convert from km^2 to m^2
% 
% stationMeanNorm = stationMean./stationArea;
% stationMaxNorm  = stationMax./stationArea;
% station100Norm  = station100./stationArea;
% 
% pntMeanNo_norm  = pntMeanNo./cellArea;
% pntMeanSP_norm  = pntMeanSP./cellArea;
% meanMaxNo_norm  = meanMaxNo./cellArea;
% meanMaxSP_norm  = meanMaxSP./cellArea;
% q100no_norm     = q100_no./cellArea;
% q100sp_norm     = q100_sp./cellArea;


%% Make preliminary plots of agreement between obs and modeled discharge
figure;
loglog(stationMean,  noMean_pnt,'*r');
%loglog(stationMeanNorm,  pntMeanNo_norm,'*r');
hold on;
loglog(stationMean, spMean_pnt, '*b');
%loglog(stationMeanNorm, pntMeanSP_norm, '*b');
hline=refline(1,0);
hold off;
hline.Color='k';
legend('noSP','SP');
htitle=title('Annual Mean Discharge');
ylabel('CaMa Discharge (m^3/s)');
xlabel('GRDC Discharge (m^3/s)');
ylim([1 10^5]);
xlim([1 10^5]);

% figure;
% loglog(stationMax,  meanMaxNo,'*r');
% hold on;
% loglog(stationMax, meanMaxSP, '*b');
% hline=refline(1,0);
% hold off;
% hline.Color='k';
% legend('noSP','SP');
% title('Annual Maximum Discharge');
% ylabel('CaMa Discharge (m^3/s)');
% xlabel('GRDC Discharge (m^3/s)');
% ylim([1 10^6]);
% xlim([1 10^6]);
% 
% figure;
% loglog(station100,  q100_no,'*r');
% hold on;
% loglog(station100, q100_sp, '*b');
% hline=refline(1,0);
% hold off;
% hline.Color='k';
% legend('noSP','SP');
% title('100-Year Discharge');
% ylabel('CaMa Discharge (m^3/s)');
% xlabel('GRDC Discharge (m^3/s)');
% ylim([1 10^6]);
% xlim([1 10^6]);

%%%%-----Make plots for h13 stations----%%%%%

% figure;
% rmean_no = corrcoef(stationMean(index),pntMeanNo(index));
% rmean_sp = corrcoef(stationMean(index),pntMeanSP(index));
% 
% loglog(stationMean(index),  pntMeanNo(index),'*r');
% hold on;
% loglog(stationMean(index), pntMeanSP(index), '*b');
% hold off;
% legend('noSP','SP');
% title('Annual Mean Discharge - H13 stations');
% ylabel('CaMa Discharge (m^3/s)');
% xlabel('GRDC Discharge (m^3/s)');
% %ylim([1 10^5]);
% %xlim([1 10^5]);
% hline=refline(1,0);
% hline.Color='k';
% % tno=text(2, 15000,['noSP Correlation Coefficient:    ', num2str(rmean_no(1,2))]);
% % tno.FontSize = 12;
% % tsp=text(2, 9000,['SP Correlation Coefficient:        ', num2str(rmean_sp(1,2))]);
% % tsp.FontSize = 12;
% 
% figure;
% rmax_no = corrcoef(stationMax(index),meanMaxNo(index));
% rmax_sp = corrcoef(stationMax(index),meanMaxSP(index));
% 
% loglog(stationMax(index),  meanMaxNo(index),'*r');
% hold on;
% loglog(stationMax(index), meanMaxSP(index), '*b');
% hold off;
% legend('noSP','SP');
% title('Annual Maximum Discharge - H13 stations');
% ylabel('CaMa Discharge (m^3/s)');
% xlabel('GRDC Discharge (m^3/s)');
% % ylim([1 10^6]);
% % xlim([1 10^6]);
% hline=refline(1,0);
% hline.Color='k';
% 
% figure;
% r100_no = corrcoef(station100(index),q100_no(index));
% r100_sp = corrcoef(station100(index),q100_sp(index));
% 
% loglog(station100(index),  q100_no(index),'*r');
% hold on;
% loglog(station100(index), q100_sp(index), '*b');
% hold off;
% legend('noSP','SP');
% title('100-Year Discharge - H13 stations');
% ylabel('CaMa Discharge (m^3/s)');
% xlabel('GRDC Discharge (m^3/s)');
% % ylim([1 10^6]);
% % xlim([1 10^6]);
% hline=refline(1,0);
% hline.Color='k';
% 

% figure;
% coast = load('coast.mat'); plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k');
% axis('xy','equal',[-180 180 -90 90]);
% hold on;
% plot(stationlon(index),stationlat(index),'.b','MarkerSize',6)
% for t=1:numel(index)
%     tname=text(stationlon(index(t)),stationlat(index(t)),stationRiver{index(t)});
%     tname.FontSize=8;
% end


