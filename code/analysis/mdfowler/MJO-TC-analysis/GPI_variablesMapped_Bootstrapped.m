% Program creates MJO phase composite maps of variables used to compute GPI
% - wind shear, relative humidity, vorticity, and potential intensity - to
% understand better the propagation of those terms, as we do OLR. 
%
% The variables are based on the bootstrapped samples used in
% GPIdecomp_barCharts_Bootstrapped.m 
%
% Meg Fowler, 2019-03-21
% Mount gp_fuse first - pritchnode.ps.uci.edu:/fast/mdfowler (for now) 

testFile  = '/Users/meganfowler/gp_fuse/TC/GPI_daily/Bootstrap_testResults/GPI_DailyGPIvars_MJOphase_1-OMI_iBoot1.nc';

lon = ncread(testFile,'lon');
lat = ncread(testFile,'lat');

load('/Users/meganfowler/gp_fuse/TC/GPI_daily/bootstrap_scripts/GPIvariables-notTerms-_bootstrapped_100.mat');
coast = load('coast.mat');


%% Average over time period and bootstraps for each variable... 

meanVort = squeeze(nanmean(absVort_boot,1));
meanVort = squeeze(nanmean(meanVort,3));

meanRH   = squeeze(nanmean(RH700_boot,1));
meanRH   = squeeze(nanmean(meanRH,3));

meanVshear = squeeze(nanmean(Vshear_boot,1));
meanVshear = squeeze(nanmean(meanVshear,3)); 

meanPI     = squeeze(nanmean(Vmax_boot,1));
meanPI     = squeeze(nanmean(meanPI,3)); 

%% Plots as made for OLR 
colormap('jet')

% ----- RH700 ----- %
figure;
c=-12:.1:12;
for iPhase=1:8
    rhPhs = squeeze(meanRH(:,:,iPhase))-squeeze(nanmean(meanRH,3));
    
    subplot(8,1,iPhase);
    contourf(lon,lat,rhPhs',c,'LineColor','none');
    title(['[Jun-Nov] Phase ',num2str(iPhase)],'FontSize',14);
    hold on;
    plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k', 'LineWidth',2);
    axis('xy','equal',[70 200, -2, 35])
    if iPhase==8
        cBar=colorbar('Ticks',[-12,-6,0,6,12],'FontSize',12);
        cBar.Label.String='700mb RH anomaly';
    end
    colormap('jet')
    caxis([min(c) max(c)]);
    
    %Add lines for each region
    line([100,100],[-2 35],'Color','white','linewidth',1.5)
    line([120,120],[-2 35],'Color','white','linewidth',1.5)
    line([130,130],[-2 35],'Color','white','linewidth',1.5)
    line([140,140],[-2 35],'Color','white','linewidth',1.5)
    line([160,160],[-2 35],'Color','white','linewidth',1.5)
    line([180,180],[-2 35],'Color','white','linewidth',1.5)
    
end
fig = gcf;
fig.PaperPositionMode = 'auto';
print('~/Documents/Irvine/TCs/Figures/phaseAnomalies_RH','-depsc');
    

% ----- Vshear ----- %
figure;
c=-6:.06:6;
for iPhase=1:8
    shearPhs = squeeze(meanVshear(:,:,iPhase))-squeeze(nanmean(meanVshear,3));
    
    subplot(8,1,iPhase);
    contourf(lon,lat,shearPhs',c,'LineColor','none');
    title(['[Jun-Nov] Phase ',num2str(iPhase)],'FontSize',14);
    hold on;
    plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k', 'LineWidth',2);
    axis('xy','equal',[70 200, -2, 35])
    if iPhase==8
        cBar=colorbar('Ticks',[-6,-3,0,3,6],'FontSize',12);
        cBar.Label.String='Wind shear anomaly';
    end
    colormap('jet')
    caxis([min(c) max(c)]);
    
    %Add lines for each region
    line([100,100],[-2 35],'Color','white','linewidth',1.5)
    line([120,120],[-2 35],'Color','white','linewidth',1.5)
    line([130,130],[-2 35],'Color','white','linewidth',1.5)
    line([140,140],[-2 35],'Color','white','linewidth',1.5)
    line([160,160],[-2 35],'Color','white','linewidth',1.5)
    line([180,180],[-2 35],'Color','white','linewidth',1.5)
end
fig = gcf;
fig.PaperPositionMode = 'auto';
print('~/Documents/Irvine/TCs/Figures/phaseAnomalies_Shear','-depsc');


% ----- Vorticity ----- %
figure;
c=-3e-6:.03e-6:3e-6;
for iPhase=1:8
    vortPhs = squeeze(meanVort(:,:,iPhase))-squeeze(nanmean(meanVort,3));
    
    subplot(8,1,iPhase);
    contourf(lon,lat,vortPhs',c,'LineColor','none');
    title(['[Jun-Nov] Phase ',num2str(iPhase)],'FontSize',14);
    hold on;
    plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k', 'LineWidth',2);
    axis('xy','equal',[70 200, -2, 35])
    if iPhase==8
        cBar=colorbar('Ticks',[-3e-6,0,3e-6],'FontSize',12);
        cBar.Label.String='Vorticity anomaly';
    end
    colormap('jet')
    caxis([min(c) max(c)]);
    
    %Add lines for each region
    line([100,100],[-2 35],'Color','white','linewidth',1.5)
    line([120,120],[-2 35],'Color','white','linewidth',1.5)
    line([130,130],[-2 35],'Color','white','linewidth',1.5)
    line([140,140],[-2 35],'Color','white','linewidth',1.5)
    line([160,160],[-2 35],'Color','white','linewidth',1.5)
    line([180,180],[-2 35],'Color','white','linewidth',1.5)
end
fig = gcf;
fig.PaperPositionMode = 'auto';
print('~/Documents/Irvine/TCs/Figures/phaseAnomalies_Vort','-depsc');


% ----- Potential Intensity ----- %
figure;
c=-4:.04:4;
for iPhase=1:8
    piPhs = squeeze(meanPI(:,:,iPhase))-squeeze(nanmean(meanPI,3));
    
    subplot(8,1,iPhase);
    contourf(lon,lat,piPhs',c,'LineColor','none');
    title(['[Jun-Nov] Phase ',num2str(iPhase)],'FontSize',14);
    hold on;
    plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k', 'LineWidth',2);
    axis('xy','equal',[70 200, -2, 35])
    if iPhase==8
        cBar=colorbar('FontSize',12);
        cBar.Label.String='Potential Intensity anomaly';
    end
    colormap('jet')
    caxis([min(c) max(c)]);
    
    %Add lines for each region
    line([100,100],[-2 35],'Color','white','linewidth',1.5)
    line([120,120],[-2 35],'Color','white','linewidth',1.5)
    line([130,130],[-2 35],'Color','white','linewidth',1.5)
    line([140,140],[-2 35],'Color','white','linewidth',1.5)
    line([160,160],[-2 35],'Color','white','linewidth',1.5)
    line([180,180],[-2 35],'Color','white','linewidth',1.5)
end
fig = gcf;
fig.PaperPositionMode = 'auto';
print('~/Documents/Irvine/TCs/Figures/phaseAnomalies_PotIntensity','-depsc');


