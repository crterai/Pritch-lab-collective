% Program creates MJO phase composite maps of variables used to compute GPI
% - wind shear, relative humidity, vorticity, and potential intensity - to
% understand better the propagation of those terms, as we do OLR. 
%
% The variables are based on the bootstrapped samples used in
% GPIdecomp_barCharts_Bootstrapped.m 
%
% Meg Fowler, 2019-03-21
% Mount gp_fuse first - sshfs mdfowler@pritchnode.ps.uci.edu:/data22/pritchard2/mdfowler/ ./gp_fuse/


%% Read in lat/lon from extra file
SSTfile  = '/Volumes/MyPassport/Data/TCs/globalSST_1983-2013_dailyAvg.nc';
% lon = ncread(SSTfile,'longitude');
% lat = ncread(SSTfile,'latitude');
SSTfile = '/Users/meganfowler/gp_fuse/TC/GPI_daily/DailyGPIvars_MJOphase_3-RMM.nc';
lon = ncread(SSTfile,'lon');
lat = ncread(SSTfile,'lat');

coast = load('coast.mat');

%% Read in GPI variables and carry out decomposition 
GPIdir   = '/Volumes/MyPassport/Data/TCs/GPI/ERA-I/';

%% Define GPI "climatology" as average of all 8 phases 

GPI_mjo     = NaN(numel(lon),numel(lat),365,8);   %Original GPI for mjo phases 
absVort_mjo = NaN(numel(lon),numel(lat),365,8);
RH700_mjo   = NaN(numel(lon),numel(lat),365,8);
Vmax_mjo    = NaN(numel(lon),numel(lat),365,8);
Vshear_mjo  = NaN(numel(lon),numel(lat),365,8);

U   = NaN(365, numel(lon),numel(lat), 2);
V   = NaN(365, numel(lon),numel(lat), 2);
SST = NaN(365, numel(lon),numel(lat));

for iPhase=1:8
    
    % -- Read in variables for GPI terms -- %
    MJOfile = [GPIdir,'GPI_DailyGPIvars_MJOphase_',num2str(iPhase),'-RMM.nc'];  %RMM index
    %MJOfile = [GPIdir,'GPI_DailyGPIvars_MJOphase',num2str(iPhase),'.nc'];   %OMI index
%     
%     absVort_mjo(:,:,:,iPhase) = ncread(MJOfile,'ABSVORT'); %[lon lat time]
%     RH700_mjo(:,:,:,iPhase)   = ncread(MJOfile,'RELHUM');
     Vmax_mjo(:,:,:,iPhase)    = ncread(MJOfile,'VMAX');
%     Vshear_mjo(:,:,:,iPhase)  = ncread(MJOfile,'VSHEAR');
%     GPI_mjo(:,:,:,iPhase) = ncread(MJOfile,'GPI'); 
% 
%     absVort_mjo(absVort_mjo<=-32767) = NaN;
%     RH700_mjo(RH700_mjo<=-32767) = NaN;
     Vmax_mjo(Vmax_mjo<=-32767) = NaN;
%     Vshear_mjo(Vshear_mjo<=-32767) = NaN;
%     GPI_mjo(GPI_mjo<=-32767) = NaN;
    
    % -- Read in variables that went into the GPI terms -- %
    varFile = ['~/gp_fuse/TC/GPI_daily/DailyGPIvars_MJOphase_',num2str(iPhase),'-RMM.nc']; %RMM index 
    %varFile = ['~/gp_fuse/TC/GPI_daily/DailyGPIvars_MJOphase',num2str(iPhase),'.nc']; %OMI index 
    
    U(:,:,:,:,iPhase) = ncread(varFile,'U'); %[m/s] 
    V(:,:,:,:,iPhase) = ncread(varFile,'V'); 
    SST(:,:,:,iPhase) = ncread(varFile,'SST'); %[K]
    
    U(U<=-32767) = NaN;
    V(V<=-32767) = NaN;
    SST(SST<=-32767) = NaN;
    
end      

levU = ncread(varFile,'levU');

%% Average over time period and bootstraps for each variable... 
iJune = 152; %Day of year that June starts
iNov  = 334; %Day of year that November ends 

% meanVort = squeeze(nanmean(absVort_mjo(:,:,iJune:iNov,:),3));
% meanRH   = squeeze(nanmean(RH700_mjo(:,:,iJune:iNov,:),3));
% meanVshear = squeeze(nanmean(Vshear_mjo(:,:,iJune:iNov,:),3));
 meanPI     = squeeze(nanmean(Vmax_mjo(:,:,iJune:iNov,:),3));

meanU   = squeeze(nanmean(U(iJune:iNov,:,:,:,:),1)); 
meanV   = squeeze(nanmean(V(iJune:iNov,:,:,:,:),1));
meanSST = squeeze(nanmean(SST(iJune:iNov,:,:,:),1));

%% Plots as made for OLR 
colormap('jet')

% ----- RH700 ----- %
% figure;
% c=-12:.1:12;
% for iPhase=1:8
%     rhPhs = squeeze(meanRH(:,:,iPhase))-squeeze(nanmean(meanRH,3));
%     
%     subplot(8,1,iPhase);
%     contourf(lon,lat,rhPhs',c,'LineColor','none');
%     title(['[Jun-Nov] Phase ',num2str(iPhase)],'FontSize',14);
%     hold on;
%     plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k', 'LineWidth',2);
%     axis('xy','equal',[70 200, -35, 35])
%     if iPhase==8
%         cBar=colorbar('Ticks',[-12,-6,0,6,12],'FontSize',12);
%         cBar.Label.String='700mb RH anomaly';
%     end
%     colormap('jet')
%     caxis([min(c) max(c)]);
%     
%     %Add lines for each region
%     line([100,100],[-35 35],'Color','white','linewidth',1.5)
%     line([120,120],[-35 35],'Color','white','linewidth',1.5)
%     line([130,130],[-35 35],'Color','white','linewidth',1.5)
%     line([140,140],[-35 35],'Color','white','linewidth',1.5)
%     line([160,160],[-35 35],'Color','white','linewidth',1.5)
%     line([180,180],[-35 35],'Color','white','linewidth',1.5)
%     
% end
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('-painters','-depsc','~/Documents/Irvine/TCs/Figures/phaseAnomalies_RH-RMM');
%     
% 
% % ----- Vshear ----- %
% figure;
% c=-6:.06:6;
% for iPhase=1:8
%     shearPhs = squeeze(meanVshear(:,:,iPhase))-squeeze(nanmean(meanVshear,3));
%     
%     subplot(8,1,iPhase);
%     contourf(lon,lat,shearPhs',c,'LineColor','none');
%     title(['[Jun-Nov] Phase ',num2str(iPhase)],'FontSize',14);
%     hold on;
%     plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k', 'LineWidth',2);
%     axis('xy','equal',[70 200, -35, 35])
%     if iPhase==8
%         cBar=colorbar('Ticks',[-6,-3,0,3,6],'FontSize',12);
%         cBar.Label.String='Wind shear anomaly';
%     end
%     colormap('jet')
%     caxis([min(c) max(c)]);
%     
%     %Add lines for each region
%     line([100,100],[-35 35],'Color','white','linewidth',1.5)
%     line([120,120],[-35 35],'Color','white','linewidth',1.5)
%     line([130,130],[-35 35],'Color','white','linewidth',1.5)
%     line([140,140],[-35 35],'Color','white','linewidth',1.5)
%     line([160,160],[-35 35],'Color','white','linewidth',1.5)
%     line([180,180],[-35 35],'Color','white','linewidth',1.5)
% end
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('~/Documents/Irvine/TCs/Figures/phaseAnomalies_Shear-RMM','-depsc');
% 
% 
% % ----- Vorticity ----- %
% figure;
% c=-3e-6:.03e-6:3e-6;
% for iPhase=1:8
%     vortPhs = squeeze(meanVort(:,:,iPhase))-squeeze(nanmean(meanVort,3));
%     
%     subplot(8,1,iPhase);
%     contourf(lon,lat,vortPhs',c,'LineColor','none');
%     title(['[Jun-Nov] Phase ',num2str(iPhase)],'FontSize',14);
%     hold on;
%     plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k', 'LineWidth',2);
%     axis('xy','equal',[70 200, -35, 35])
%     if iPhase==8
%         cBar=colorbar('Ticks',[-3e-6,0,3e-6],'FontSize',12);
%         cBar.Label.String='Vorticity anomaly';
%     end
%     colormap('jet')
%     caxis([min(c) max(c)]);
%     
%     %Add lines for each region
%     line([100,100],[-35 35],'Color','white','linewidth',1.5)
%     line([120,120],[-35 35],'Color','white','linewidth',1.5)
%     line([130,130],[-35 35],'Color','white','linewidth',1.5)
%     line([140,140],[-35 35],'Color','white','linewidth',1.5)
%     line([160,160],[-35 35],'Color','white','linewidth',1.5)
%     line([180,180],[-35 35],'Color','white','linewidth',1.5)
% end
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('~/Documents/Irvine/TCs/Figures/phaseAnomalies_Vort-RMM','-depsc');
% 

% ----- Potential Intensity ----- %
figure;
c=-6:.06:6;
for iPhase=1:8
    piPhs = squeeze(meanPI(:,:,iPhase))-squeeze(nanmean(meanPI,3));
    
    subplot(8,1,iPhase);
    contourf(lon,lat,piPhs',c,'LineColor','none');
    title(['[Jun-Nov] Phase ',num2str(iPhase)],'FontSize',14);
    hold on;
    plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k', 'LineWidth',2);
    axis('xy','equal',[70 200, -35, 35])
    if iPhase==8
        cBar=colorbar('FontSize',12);
        cBar.Label.String='Potential Intensity anomaly';
    end
    colormap('jet')
    caxis([min(c) max(c)]);
    
    %Add lines for each region
    line([100,100],[-35 35],'Color','white','linewidth',1.5)
    line([120,120],[-35 35],'Color','white','linewidth',1.5)
    line([130,130],[-35 35],'Color','white','linewidth',1.5)
    line([140,140],[-35 35],'Color','white','linewidth',1.5)
    line([160,160],[-35 35],'Color','white','linewidth',1.5)
    line([180,180],[-35 35],'Color','white','linewidth',1.5)
end
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-painters','-depsc','~/Documents/Irvine/TCs/Figures/phaseAnomalies_PotIntensity-RMM');


% ----- SST ----- %
figure;
c=-0.5:.005:0.5;
for iPhase=1:8
    sstPhs = squeeze(meanSST(:,:,iPhase))-squeeze(nanmean(meanSST,3));
    
    subplot(8,1,iPhase);
    contourf(lon,lat,sstPhs',c,'LineColor','none');
    title(['[Jun-Nov] Phase ',num2str(iPhase)],'FontSize',14);
    hold on;
    plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k', 'LineWidth',2);
    axis('xy','equal',[70 200, -35, 35])
    if iPhase==8
        cBar=colorbar('FontSize',12);
        cBar.Label.String='SST anomaly';
    end
    colormap('jet')
    caxis([min(c) max(c)]);
    
    %Add lines for each region
    line([100,100],[-35 35],'Color','white','linewidth',1.5)
    line([120,120],[-35 35],'Color','white','linewidth',1.5)
    line([130,130],[-35 35],'Color','white','linewidth',1.5)
    line([140,140],[-35 35],'Color','white','linewidth',1.5)
    line([160,160],[-35 35],'Color','white','linewidth',1.5)
    line([180,180],[-35 35],'Color','white','linewidth',1.5)
end
fig = gcf;
fig.PaperPositionMode = 'auto';
print('-painters','-depsc','~/Documents/Irvine/TCs/Figures/phaseAnomalies_SST-RMM');


% % ----- U250 ----- %
% figure;
% c=-6:.06:6;
% for iPhase=1:8
%     u250Phs = squeeze(meanU(:,:,1,iPhase))-squeeze(nanmean(meanU(:,:,1,:),4));
%     
%     subplot(8,1,iPhase);
%     contourf(lon,lat,u250Phs',c,'LineColor','none');
%     title(['[Jun-Nov] Phase ',num2str(iPhase)],'FontSize',14);
%     hold on;
%     plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k', 'LineWidth',2);
%     axis('xy','equal',[70 200, -35, 35])
%     if iPhase==8
%         cBar=colorbar('FontSize',12);
%         cBar.Label.String='U250 anomaly';
%     end
%     colormap('jet')
%     caxis([min(c) max(c)]);
%     
%     %Add lines for each region
%     line([100,100],[-35 35],'Color','white','linewidth',1.5)
%     line([120,120],[-35 35],'Color','white','linewidth',1.5)
%     line([130,130],[-35 35],'Color','white','linewidth',1.5)
%     line([140,140],[-35 35],'Color','white','linewidth',1.5)
%     line([160,160],[-35 35],'Color','white','linewidth',1.5)
%     line([180,180],[-35 35],'Color','white','linewidth',1.5)
% end
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('-painters','-depsc','~/Documents/Irvine/TCs/Figures/phaseAnomalies_U250-RMM');
% 
% 
% % ----- U850 ----- %
% figure;
% c=-5:.05:5;
% for iPhase=1:8
%     u250Phs = squeeze(meanU(:,:,2,iPhase))-squeeze(nanmean(meanU(:,:,2,:),4));
%     
%     subplot(8,1,iPhase);
%     contourf(lon,lat,u250Phs',c,'LineColor','none');
%     title(['[Jun-Nov] Phase ',num2str(iPhase)],'FontSize',14);
%     hold on;
%     plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k', 'LineWidth',2);
%     axis('xy','equal',[70 200, -35, 35])
%     if iPhase==8
%         cBar=colorbar('FontSize',12);
%         cBar.Label.String='U850 anomaly';
%     end
%     colormap('jet')
%     caxis([min(c) max(c)]);
%     
%     %Add lines for each region
%     line([100,100],[-35 35],'Color','white','linewidth',1.5)
%     line([120,120],[-35 35],'Color','white','linewidth',1.5)
%     line([130,130],[-35 35],'Color','white','linewidth',1.5)
%     line([140,140],[-35 35],'Color','white','linewidth',1.5)
%     line([160,160],[-35 35],'Color','white','linewidth',1.5)
%     line([180,180],[-35 35],'Color','white','linewidth',1.5)
% end
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print('-painters','-depsc','~/Documents/Irvine/TCs/Figures/phaseAnomalies_U850-RMM');


%% Maps of background states  (Phase 1-8 average) 

% --- U850 --- %
figure; 
contourf(lon,lat,squeeze(nanmean(meanU(:,:,2,:),4))',-10:.1:10,'LineColor','none');
title('[Jun-Nov] Average U850 (Phase1-8 Mean)','FontSize',14);
hold on;
plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k', 'LineWidth',2);
axis('xy','equal',[70 200, -35, 35])
cBar=colorbar('FontSize',12);
cBar.Label.String='U850 [m/s]';
colormap('jet')
caxis([-10 10]);
%Add lines for each region
line([100,100],[-35 35],'Color','white','linewidth',1.5)
line([120,120],[-35 35],'Color','white','linewidth',1.5)
line([130,130],[-35 35],'Color','white','linewidth',1.5)
line([140,140],[-35 35],'Color','white','linewidth',1.5)
line([160,160],[-35 35],'Color','white','linewidth',1.5)
line([180,180],[-35 35],'Color','white','linewidth',1.5)
%Save Figure 
fig = gcf;
fig.PaperPositionMode = 'auto';
print('~/Documents/Irvine/TCs/Figures/phaseMean_U850','-depsc');

% --- U250 --- %
figure; 
contourf(lon,lat,squeeze(nanmean(meanU(:,:,1,:),4))',-25:.25:25,'LineColor','none');
title('[Jun-Nov] Average U250 (Phase1-8 Mean)','FontSize',14);
hold on;
plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k', 'LineWidth',2);
axis('xy','equal',[70 200, -35, 35])
cBar=colorbar('FontSize',12);
cBar.Label.String='U250 [m/s]';
colormap('jet')
caxis([-25 25]);
%Add lines for each region
line([100,100],[-35 35],'Color','white','linewidth',1.5)
line([120,120],[-35 35],'Color','white','linewidth',1.5)
line([130,130],[-35 35],'Color','white','linewidth',1.5)
line([140,140],[-35 35],'Color','white','linewidth',1.5)
line([160,160],[-35 35],'Color','white','linewidth',1.5)
line([180,180],[-35 35],'Color','white','linewidth',1.5)
%Save Figure 
fig = gcf;
fig.PaperPositionMode = 'auto';
print('~/Documents/Irvine/TCs/Figures/phaseMean_U250','-depsc');

% --- SST --- %
figure; 
contourf(lon,lat,squeeze(nanmean(meanSST,3))',290:0.2:305,'LineColor','none');
title('[Jun-Nov] Average SST (Phase1-8 Mean)','FontSize',14);
hold on;
plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k', 'LineWidth',2);
axis('xy','equal',[70 200, -35, 35])
cBar=colorbar('FontSize',12);
cBar.Label.String='SST [K]';
colormap('jet')
caxis([290 305]);
%Add lines for each region
line([100,100],[-35 35],'Color','white','linewidth',1.5)
line([120,120],[-35 35],'Color','white','linewidth',1.5)
line([130,130],[-35 35],'Color','white','linewidth',1.5)
line([140,140],[-35 35],'Color','white','linewidth',1.5)
line([160,160],[-35 35],'Color','white','linewidth',1.5)
line([180,180],[-35 35],'Color','white','linewidth',1.5)
%Save Figure 
fig = gcf;
fig.PaperPositionMode = 'auto';
print('~/Documents/Irvine/TCs/Figures/phaseMean_SST','-depsc');


% --- Wind Shear --- %
figure; 
contourf(lon,lat,squeeze(nanmean(meanVshear,3))',0:0.4:40,'LineColor','none');
title('[Jun-Nov] Average Wind Shear (Phase1-8 Mean)','FontSize',14);
hold on;
plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k', 'LineWidth',2);
axis('xy','equal',[70 200, -35, 35])
cBar=colorbar('FontSize',12);
cBar.Label.String='Wind Shear [m/s]';
colormap('jet')
caxis([0 40]);
%Add lines for each region
line([100,100],[-35 35],'Color','white','linewidth',1.5)
line([120,120],[-35 35],'Color','white','linewidth',1.5)
line([130,130],[-35 35],'Color','white','linewidth',1.5)
line([140,140],[-35 35],'Color','white','linewidth',1.5)
line([160,160],[-35 35],'Color','white','linewidth',1.5)
line([180,180],[-35 35],'Color','white','linewidth',1.5)
%Save Figure 
fig = gcf;
fig.PaperPositionMode = 'auto';
print('~/Documents/Irvine/TCs/Figures/phaseMean_Vshear','-depsc');

%% Try regional bar charts, as done for GPI 

%Get regional averages of each term 
lonsLeft  = [100,120,130,140,160];
lonsRight = [120,130,140,160,180];
ilats = find(lat>=5 & lat<=20);

% for iPhase=1:8
%     dU850_pct(:,:,iPhase)   = (squeeze(meanU(:,:,2,iPhase))-squeeze(nanmean(meanU(:,:,2,:),4)))./squeeze(nanmean(meanU(:,:,2,:),4)); 
%     dU250_pct(:,:,iPhase)   = (squeeze(meanU(:,:,1,iPhase))-squeeze(nanmean(meanU(:,:,1,:),4)))./squeeze(nanmean(meanU(:,:,1,:),4)); 
%     dSST_pct(:,:,iPhase)    = (squeeze(meanSST(:,:,iPhase))-squeeze(nanmean(meanSST,3)))./squeeze(nanmean(meanSST,3)); 
%     
%     for iReg=1:5
%        ilons = find(lon>=lonsLeft(iReg) & lon<lonsRight(iReg)); 
% 
%        dU850(iPhase,iReg)     = nanmean(nanmean(dU850_pct(ilons,ilats,iPhase)))*100;
%        dU250(iPhase,iReg)     = nanmean(nanmean(dU250_pct(ilons,ilats,iPhase)))*100;
%        dSST(iPhase,iReg)      = nanmean(nanmean(dSST_pct(ilons,ilats,iPhase)))*100;
% 
%        clearvars ilons 
%     end
% end
% 
% 
% figure;
% c            = categorical({'Phase 1', 'Phase 2', 'Phase 3', 'Phase 4', 'Phase 5', 'Phase 6', 'Phase 7', 'Phase 8'});
% regionLabels = {'100?-120?','120?-130?','130?-140?','140?-160?','160?-180?'};

% for iGPI=1:5   
%     phsRegions = [dU850(1,iGPI),dU250(1,iGPI),dSST(1,iGPI);...
%                   dU850(2,iGPI),dU250(2,iGPI),dSST(2,iGPI);...
%                   dU850(3,iGPI),dU250(3,iGPI),dSST(3,iGPI);...
%                   dU850(4,iGPI),dU250(4,iGPI),dSST(4,iGPI);...
%                   dU850(5,iGPI),dU250(5,iGPI),dSST(5,iGPI);...
%                   dU850(6,iGPI),dU250(6,iGPI),dSST(6,iGPI);...
%                   dU850(7,iGPI),dU250(7,iGPI),dSST(7,iGPI);...
%                   dU850(8,iGPI),dU250(8,iGPI),dSST(8,iGPI)];
%         
%     %Only add x-axis regions if bottom plot
%     subplot(3,2,iGPI)     
%     bar(phsRegions);
%     %ylim([-1.5, 1.5]);
%     xlabel('MJO Phase [OMI]','FontSize',16)
%     title(['Region: ',regionLabels{iGPI}],'fontsize',16);
%     ax = gca;
%     ax.FontSize = 14;
%     
%     
%     %Only add legend to top plot
%     if iGPI==5
%         legend({'dU850','dU250','dSST'},'location','southeast','fontsize',12);
%     end
%     
%        
%     clearvars phsRegions
%     %print(['~/Desktop/TestBarChart_Phs',num2str(iGPI),'.pdf'],'-dpdf','-bestfit')
% end


% --- What if I just plot the mean wind for each phase, instead of 'anomaly' --- %

% Apply ocean mask for making bar charts 
mask = double(isfinite(meanSST));
mask(mask==0) = NaN;
U850_mask = squeeze(meanU(:,:,2,:)).*mask;
U250_mask = squeeze(meanU(:,:,1,:)).*mask;


for iPhase=1:8
    for iReg=1:5
       ilons = find(lon>=lonsLeft(iReg) & lon<lonsRight(iReg)); 

       U850_reg(iPhase,iReg)     = nanmean(nanmean(squeeze(U850_mask(ilons,ilats,iPhase))));
       U250_reg(iPhase,iReg)     = nanmean(nanmean(squeeze(U250_mask(ilons,ilats,iPhase))));
       SST_reg(iPhase,iReg)      = nanmean(nanmean(squeeze(meanSST(ilons,ilats,iPhase))));

       clearvars ilons 
    end
end


figure;
c            = categorical({'Phase 1', 'Phase 2', 'Phase 3', 'Phase 4', 'Phase 5', 'Phase 6', 'Phase 7', 'Phase 8'});
regionLabels = {'100?-120?','120?-130?','130?-140?','140?-160?','160?-180?'};

for iGPI=1:5   
    phsRegions = [U850_reg(1,iGPI), U250_reg(1,iGPI);...
                  U850_reg(2,iGPI), U250_reg(2,iGPI);...
                  U850_reg(3,iGPI), U250_reg(3,iGPI);...
                  U850_reg(4,iGPI), U250_reg(4,iGPI);...
                  U850_reg(5,iGPI), U250_reg(5,iGPI);...
                  U850_reg(6,iGPI), U250_reg(6,iGPI);...
                  U850_reg(7,iGPI), U250_reg(7,iGPI);...
                  U850_reg(8,iGPI), U250_reg(8,iGPI)];
    
    mean850 = nanmean(U850_reg(:,iGPI),1);
    mean250 = nanmean(U250_reg(:,iGPI),1);
              
    %Only add x-axis regions if bottom plot
    subplot(3,2,iGPI)     
    b = bar(phsRegions);
    %ylim([-1.5, 1.5]);
    xlabel('MJO Phase [OMI]','FontSize',16)
    ylabel('Wind Speed [m/s]');
    title(['Region: ',regionLabels{iGPI}],'fontsize',16);
    ax = gca;
    ax.FontSize = 14;
    
    %Control color of bars
    b(1).FaceColor = [0 0.4470 0.7410];
    b(2).FaceColor = [0.4660 0.6740 0.1880];
    
    %Add horizontal line for phase mean wind speeds in this region 
    hold on; 
    plot(xlim,[mean850 mean850],'b--','LineWidth',2)
    plot(xlim,[mean250 mean250],'--','Color',[0, (204/255),0],'LineWidth',2)
    hold off;
    
    %Only add legend to top plot
    if iGPI==5
        legend({'U850','U250'},'location','southeast','fontsize',12);
    end
    
       
    clearvars phsRegions
    %print(['~/Desktop/TestBarChart_Phs',num2str(iGPI),'.pdf'],'-dpdf','-bestfit')
end



