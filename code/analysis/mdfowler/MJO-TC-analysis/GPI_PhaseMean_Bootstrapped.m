% Map Phase 1-8 mean GPI 
%
% 24 July 2019
% Meg Fowler
% 
% Mount gp_fuse first - pritchnode.ps.uci.edu:/fast/mdfowler (for now) 

%% Read in data

GPI_dir = '/Users/meganfowler/gp_fuse/TC/GPI_daily/Bootstrap_testResults/';

%Define lat and lon 
testFile = [GPI_dir,'GPI_DailyGPIvars_MJOphase_',num2str(1),'-OMI_iBoot',num2str(1),'.nc'];
lat      = ncread(testFile,'lat');
lon      = ncread(testFile,'lon');

%Define empty arrays 
allGPI = nan(8,100,length(lon),length(lat),365);

%Read in  GPI for all phases and bootstraps 
for iPhase = 1:8
    
    for iBoot = 1:99
        fileName = [GPI_dir,'GPI_DailyGPIvars_MJOphase_',num2str(iPhase),'-OMI_iBoot',num2str(iBoot),'.nc'];
        allGPI(iPhase,iBoot,:,:,:) = ncread(fileName,'GPI');


    end
    fprintf('Done with phase %i \n', iPhase)
end

%% Take averages 

%Average over bootstraps 
GPI_bootMean = squeeze(nanmean(allGPI,2));

%Average over TC season 
iJune = 152; %Day of year that June starts
iNov  = 334; %Day of year that November ends 

GPI_timeMean = squeeze(nanmean(GPI_bootMean(:,:,:,iJune:iNov),4));

%Average over 8 phases
GPI_phsMean = squeeze(nanmean(GPI_timeMean,1));

%% Or, take monthly averages: 


%Average over each month 
GPI_monthly(:,:,:,1) = squeeze(nanmean(GPI_bootMean(:,:,:,1:31),4));     %Jan - 31 days
GPI_monthly(:,:,:,2) = squeeze(nanmean(GPI_bootMean(:,:,:,32:59),4));    %Feb - 28 days
GPI_monthly(:,:,:,3) = squeeze(nanmean(GPI_bootMean(:,:,:,60:90),4));    %Mar - 31 days
GPI_monthly(:,:,:,4) = squeeze(nanmean(GPI_bootMean(:,:,:,91:120),4));   %Apr - 30 days
GPI_monthly(:,:,:,5) = squeeze(nanmean(GPI_bootMean(:,:,:,121:151),4));  %May - 31 days
GPI_monthly(:,:,:,6) = squeeze(nanmean(GPI_bootMean(:,:,:,152:181),4));  %Jun - 30 days
GPI_monthly(:,:,:,7) = squeeze(nanmean(GPI_bootMean(:,:,:,182:212),4));  %Jul - 31 days
GPI_monthly(:,:,:,8) = squeeze(nanmean(GPI_bootMean(:,:,:,213:243),4));  %Aug - 31 days
GPI_monthly(:,:,:,9) = squeeze(nanmean(GPI_bootMean(:,:,:,244:273),4));  %Sep - 30 days
GPI_monthly(:,:,:,10) = squeeze(nanmean(GPI_bootMean(:,:,:,274:304),4)); %Oct - 31 days
GPI_monthly(:,:,:,11) = squeeze(nanmean(GPI_bootMean(:,:,:,305:334),4)); %Nov - 30 days
GPI_monthly(:,:,:,12) = squeeze(nanmean(GPI_bootMean(:,:,:,335:365),4)); %Dec - 31 days

%Average over all phases 
GPI_monthly_phsMean = squeeze(nanmean(GPI_monthly,1));



%% Plotting; follow format for MIT model gen density plots 
% 
% minLon = find(lon==100);
% maxLon = find(lon==200);
% minLat = find(lat==30);
% maxLat = find(lat==0);
% 
% m_proj(projection,'long',[100 200],'lat',[0 30]);
% cmax=(max(max(GPI_phsMean)));
% cint=cmax/15;
% %[C1,h1]=m_contour(x(kmin:kmax),y(imin:imax),z(imin:imax,kmin:kmax),cint:cint:cmax);
% [C1,h1]=m_contour(lon(minLon:maxLon),lat(minLat:maxLat),GPI_phsMean(minLon:maxLon, minLat:maxLat)',cint:cint:cmax,'LineColor','none');   %Modified to have same colorbar
% colormap(wmap)
% if strcmp(wfill,'y')
%     hold on
%     h=m_pcolor(lon(minLon:maxLon),lat(minLat:maxLat),GPI_phsMean(minLon:maxLon, minLat:maxLat)');
%     a=colormap(wmap);
%     for i=1:8
%         a(i,:)=oceancolor;
%     end
%     colormap(a)
%     set(h,'EdgeAlpha',0,'FaceAlpha',wtrans);
%     set(h,'FaceColor','interp');
%     hold off
% end
% cb=colorbar;
% caxis([cint cmax])
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
% 
% title('Phase 1-8 Mean GPI','fontsize',12,'fontweight','bold','interpreter','none')



