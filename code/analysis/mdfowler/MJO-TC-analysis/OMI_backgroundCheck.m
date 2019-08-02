% Check SST anomalies in OMI phases - did something go wrong?
%
% Meg D. Fowler, 2017-10-26
%

%% Load daily SST and compute average over full period 
rootPath = '/Users/meganfowler/gp_fuse/TC/MJO_obs/ECMWF_ERA-Interim/2.5x2.5/'; %Path to Greenplanet ECMWF obs
SSTfile  = [rootPath, 'dailyAvg/globalSST_1983-2013_dailyAvg.nc'];

%SSTfull = ncread(SSTfile,'sst');
lon     = ncread(SSTfile,'longitude');
lat     = ncread(SSTfile,'latitude');
% SSTmean = squeeze(nanmean(SSTfull,3));  %Average over all the days 
% clearvars SSTfull SSTfile rootPath 
% 
% %% Load phase-separated SST from files 
% rootPath = '/Users/meganfowler/gp_fuse/TC/MJO_obs/ECMWF_ERA-Interim/phsPair-ForKerry/OMI/';
% fileBase = 'OMI_annual_MJOphase';
% coast = load('coast.mat');
% 
% figure;
% SSTanom = NaN(144,73,8);
% for iPhase = 1:8
%     filename = [rootPath,fileBase,num2str(iPhase),'.nc'];
%     SSTphs = ncread(filename,'SST');
%     SSTmeanPhs = squeeze(nanmean(SSTphs,1));   %Average over all the months 
%     
%     SSTanom(:,:,iPhase) = SSTmeanPhs - SSTmean; 
%     
%     subplot(8,1,iPhase);
%     contourf(lon,lat,squeeze(SSTanom(:,:,iPhase))','LineColor','none');
%     title(sprintf('SST anomaly: Phase %d [OMI]',iPhase));
%     hold on; 
%     plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k');
%     axis('xy','equal',[0 280 -40 40]);
%     colorbar;
%     caxis([-0.3 0.3])    
%     
% end

%% Determine SST anomaly relative to phase 1-8 mean... 
OMIdir = '/Users/meganfowler/gp_fuse/TC/MJO_obs/ECMWF_ERA-Interim/phsPair-ForKerry/OMI/OMI_annual_MJOphase';

for iPhs=1:8
   fileIn = [OMIdir,num2str(iPhs),'.nc']; 
   
   SST(:,:,:,iPhs) = ncread(fileIn,'SST'); %Read in monthly average SST 
end

%Compute SST anomaly from phase 1-8 average
SSTmean = squeeze(nanmean(SST,4));

figure;
coast = load('coast.mat');
cSST = -0.5:0.005:0.5;

for iPhs = 1:8 
   SSTanom = SST(:,:,:,iPhs) - SSTmean;  
   
   %Plot shear anomaly 
   subplot(8,1,iPhs);
   contourf(lon,lat,squeeze(nanmean(SSTanom(6:11,:,:),1))',cSST,'LineColor','none');
   title(sprintf('June-Nov SST Anomaly: Phase %d [OMI]',iPhs));
   hold on; 
   plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k');
   axis('xy','equal',[0 200 -40 40]);
   colorbar;
   caxis([min(cSST) max(cSST)])  
   
end

%% Create similar plots for 850 and 250 mb winds 
OMIdir = '/Users/meganfowler/gp_fuse/TC/MJO_obs/ECMWF_ERA-Interim/phsPair-ForKerry/OMI/OMI_annual_MJOphase';

for iPhs=1:8
   fileIn = [OMIdir,num2str(iPhs),'.nc']; 
   
   U(:,:,:,:,iPhs) = ncread(fileIn,'U'); %Read in daily U 
   V(:,:,:,:,iPhs) = ncread(fileIn,'V'); %Read in daily V 
end

levU = ncread(fileIn,'levU');

%Define individual levels 
U250 = squeeze(U(:,:,:,1,:));
U850 = squeeze(U(:,:,:,2,:));
V250 = squeeze(V(:,:,:,1,:));
V850 = squeeze(V(:,:,:,2,:));

%Take mean over all 8 phases 
U250mean = squeeze(nanmean(U250,4));
U850mean = squeeze(nanmean(U850,4));
V250mean = squeeze(nanmean(V250,4));
V850mean = squeeze(nanmean(V850,4));

figure(1)
figure(2)
figure(3)
figure(4) 
cU = -5:0.05:5;
for iPhs=1:8
    U250anom = U250(:,:,:,iPhs) - U250mean;  
    U850anom = U850(:,:,:,iPhs) - U850mean;  
    V250anom = V250(:,:,:,iPhs) - V250mean;  
    V850anom = V850(:,:,:,iPhs) - V850mean;  
   
   figure(1)
   %Plot U250 anomaly  
   subplot(8,1,iPhs);
   contourf(lon,lat,squeeze(nanmean(U250anom(152:334,:,:),1))',cU,'LineColor','none');
   title(sprintf('June-Nov U250: Phase %d [OMI]',iPhs));
   hold on; 
   plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k');
    axis([80 200 0 40]);
   colorbar;
   caxis([min(cU) max(cU)])  

   figure(2)
   %Plot U850 anomaly  
   subplot(8,1,iPhs);
   contourf(lon,lat,squeeze(nanmean(U850anom(152:334,:,:),1))',cU,'LineColor','none');
   title(sprintf('June-Nov U850: Phase %d [OMI]',iPhs));
   hold on; 
   plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k');
    axis([80 200 0 40]);
   colorbar;
   caxis([min(cU) max(cU)])  
   
   figure(3)
   %Plot V850 anomaly  
   subplot(8,1,iPhs);
   contourf(lon,lat,squeeze(nanmean(V850anom(152:334,:,:),1))',cU,'LineColor','none');
   title(sprintf('June-Nov V850: Phase %d [OMI]',iPhs));
   hold on; 
   plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k');
    axis([80 200 0 40]);
   colorbar;
   caxis([min(cU) max(cU)])  
  
   figure(4)
   %Plot V250 anomaly  
   subplot(8,1,iPhs);
   contourf(lon,lat,squeeze(nanmean(V250anom(152:334,:,:),1))',cU,'LineColor','none');
   title(sprintf('June-Nov V250: Phase %d [OMI]',iPhs));
   hold on; 
   plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k');
    axis([80 200 0 40]);
   colorbar;
   caxis([min(cU) max(cU)]) 
    
    
end


%% Create similar plot for vertical wind shear 
%OMIdir = '/Users/meganfowler/gp_fuse/TC/MJO_obs/ECMWF_ERA-Interim/phsPair-ForKerry/OMI/'; 
GPIdir = '/Volumes/MyPassport/Data/TCs/GPI/ERA-I/';

%Read in files and data 
for iPhs = 1:8
   fileIn = [GPIdir,'GPI_DailyGPIvars_MJOphase',num2str(iPhs),'.nc'];
   
   shear(:,:,:,iPhs) = ncread(fileIn,'VSHEAR');   %Read in daily wind shear   
end

%Compute wind shear anomaly from phase 1-8 average 
shearMean = squeeze(nanmean(shear,4));

figure;
coast = load('coast.mat');
cShear = -7:0.07:7;

for iPhs = 1:8 
   shearAnom = shear(:,:,:,iPhs) - shearMean;  
   
   %Plot shear anomaly 
   subplot(8,1,iPhs);
   contourf(lon,lat,squeeze(nanmean(shearAnom(:,:,152:334),3))',cShear,'LineColor','none');
   title(sprintf('June-Nov Wind Shear Anomaly: Phase %d [OMI]',iPhs));
   hold on; 
   plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k');
   axis([60 120 0 40]);
   colorbar;
   caxis([min(cShear) max(cShear)])  
   
end


%% Determine RH anomaly relative to Phase 1-8 mean
GPIdir = '/Volumes/MyPassport/Data/TCs/GPI/ERA-I/';

for iPhs=1:8
   fileIn = [GPIdir,'GPI_DailyGPIvars_MJOphase',num2str(iPhs),'.nc'];
   
   RH(:,:,:,iPhs) = ncread(fileIn,'RELHUM'); %Read in daily 700 mb RH 
   
end
RH(RH>=9.96e36) = NaN;

%Compute SST anomaly from phase 1-8 average
RHmean = squeeze(nanmean(RH,4));

figure;
coast = load('coast.mat');
cRH = -13:0.1:13;

for iPhs = 1:8 
   RHanom = RH(:,:,:,iPhs) - RHmean;  
   
   %Plot shear anomaly 
   subplot(8,1,iPhs);
   contourf(lon,lat,squeeze(nanmean(RHanom(:,:,152:334),3))',cRH,'LineColor','none');
   title(sprintf('June-Nov 700mb RH Anomaly: Phase %d [OMI]',iPhs));
   hold on; 
   plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k');
   axis('xy','equal',[0 200 -40 40]);
   colorbar;
   caxis([min(cRH) max(cRH)])  
   
end

%% Determine absolute vorticity anomaly relative to Phase 1-8 mean
GPIdir = '/Volumes/MyPassport/Data/TCs/GPI/ERA-I/';

for iPhs=1:8
   fileIn = [GPIdir,'GPI_DailyGPIvars_MJOphase',num2str(iPhs),'.nc'];
   
   Vort(:,:,:,iPhs) = ncread(fileIn,'ABSVORT'); %Read in daily absolute vorticity
   
end
Vort(Vort==-32767) = NaN;

%Compute SST anomaly from phase 1-8 average
Vortmean = squeeze(nanmean(Vort,4));

figure;
coast = load('coast.mat');
cV = -10e-6:10e-8:10e-6;

for iPhs = 1:8 
   VortAnom = Vort(:,:,:,iPhs) - Vortmean;  
   
   %Plot shear anomaly 
   subplot(8,1,iPhs);
   contourf(lon,lat,squeeze(nanmean(VortAnom(:,:,152:334),3))',cV,'LineColor','none');
   title(sprintf('June-Nov Abs Vorticity Anomaly: Phase %d [OMI]',iPhs));
   hold on; 
   plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k');
   axis('xy','equal',[0 200 -40 40]);
   colorbar;
   caxis([min(cV) max(cV)])  
   
end