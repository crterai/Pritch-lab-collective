% Make plot of bootstrapped background wind characteristics (bar charts
% with lines for the mean) 
%
% Meg Fowler, 2019-05-08

SSTfile = '/Users/meganfowler/gp_fuse/TC/GPI_daily/DailyGPIvars_MJOphase_3-RMM.nc';
lon  = ncread(SSTfile,'lon');
lat  = ncread(SSTfile,'lat');
levU = ncread(SSTfile,'levU');
SST  = ncread(SSTfile,'SST');
SST  = squeeze(SST(1,:,:));

%Domain in bootstrapped winds file 
ilon = find(lon>=60 & lon<=200);
ilat = find(lat>-5 & lat<40);

lon = lon(ilon);
lat = lat(ilat); 
SST = SST(ilon,ilat); 

%% Read in GPI variables for each bootstrap -- done on Greenplanet

% varDir = '/Users/meganfowler/gp_fuse/TC/GPI_daily/Bootstrap_DailyGPIvariables/';
% U = NaN(365,57,17,2,8,100);
% 
% for iBoot = 1:100
%     for iPhs = 1:8
%         tic
%         fileName = [varDir,'DailyGPIvars_MJOphase_',num2str(iPhs),'-OMI_iBoot',num2str(iBoot),'.nc'];
%         
%         U(:,:,:,:,iPhs,iBoot) = ncread(fileName,'U');
%         toc
%     end
%     fprintf('Done with bootstrap $i \n',iBoot);
% end
% U250 = squeeze(U(:,:,:,1,:,:));
% U850 = squeeze(U(:,:,:,2,:,:));
% 
% U250avg = squeeze(nanmean(U250,5));
% U850avg = squeeze(nanmean(U850,5));
% 
% U250_pctile = prctile(U250, [25,75], 5);
% U850_pctile = prctile(U850, [25,75], 5);
% 
% save('bootstrappedWinds.mat','U250avg', 'U850avg','U250_pctile','U850_pctile');

%% Ran that on Pritchnode - just read in final .mat file with winds at both levels 
load('/Users/meganfowler/gp_fuse/TC/GPI_daily/bootstrap_scripts/bootstrappedWinds.mat');

%% Get regional averages of each term 
lonsLeft  = [100,120,130,140,160];
lonsRight = [120,130,140,160,180];
ilats = find(lat>=5 & lat<=20);

iJune = 152; %Day of year that June starts
iNov  = 334; %Day of year that November ends 

%Average over select season
U250avg = squeeze(nanmean(U250avg(iJune:iNov,:,:,:),1)); 
U850avg = squeeze(nanmean(U850avg(iJune:iNov,:,:,:),1)); 
U250_pctile = squeeze(nanmean(U250_pctile(iJune:iNov,:,:,:,:),1)); 
U850_pctile = squeeze(nanmean(U850_pctile(iJune:iNov,:,:,:,:),1)); 

% Apply ocean mask  
mask = double(isfinite(SST));
mask(mask==0) = NaN;
for iPhs=1:8 
    mask_phs(:,:,iPhs) = mask; %Expand to fill all phases 
end

U850_mask = U850avg.*mask_phs;
U250_mask = U250avg.*mask_phs;
U850_25   = squeeze(U850_pctile(:,:,:,1)).*mask_phs;
U850_75   = squeeze(U850_pctile(:,:,:,2)).*mask_phs; 
U250_25   = squeeze(U250_pctile(:,:,:,1)).*mask_phs;
U250_75   = squeeze(U250_pctile(:,:,:,2)).*mask_phs; 

% Regional means 
for iPhase=1:8
    for iReg=1:5
       ilons = find(lon>=lonsLeft(iReg) & lon<lonsRight(iReg)); 

       U850_reg(iPhase,iReg)     = nanmean(nanmean(squeeze(U850_mask(ilons,ilats,iPhase))));
       U250_reg(iPhase,iReg)     = nanmean(nanmean(squeeze(U250_mask(ilons,ilats,iPhase))));
       U850_25_reg(iPhase,iReg)  = nanmean(nanmean(squeeze(U850_25(ilons,ilats,iPhase))));
       U850_75_reg(iPhase,iReg)  = nanmean(nanmean(squeeze(U850_75(ilons,ilats,iPhase))));
       U250_25_reg(iPhase,iReg)  = nanmean(nanmean(squeeze(U250_25(ilons,ilats,iPhase))));
       U250_75_reg(iPhase,iReg)  = nanmean(nanmean(squeeze(U250_75(ilons,ilats,iPhase))));

       clearvars ilons 
        
        
    end
end

           

%% Make wind bar chart plot with error bars now 
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

fig = gcf;
fig.PaperPositionMode = 'auto';
print('/Users/meganfowler/Documents/Irvine/TCs/Figures/WindBarChart_BasedOnBootstraps','-depsc');



