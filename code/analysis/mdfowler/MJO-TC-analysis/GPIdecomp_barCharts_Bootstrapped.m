% Create bar chart for GPI term contributions to full GPI, as plotted in
% map form on the AMS poster and most of my talks. 
%
% 2019-03-05 

%% Read in lat/lon from extra file
% testFile  = '/Users/meganfowler/gp_fuse/TC/GPI_daily/Bootstrap_testResults/GPI_DailyGPIvars_MJOphase_1-OMI_iBoot1.nc';
% 
% lon = ncread(testFile,'lon');
% lat = ncread(testFile,'lat');
% 
% iJune = 152; %Day of year that June starts
% iNov  = 334; %Day of year that November ends 
% % 
% % 
% % %% Read in GPI variables and carry out decomposition 
% % %GPIdir   = '/Volumes/MyPassport/Data/TCs/GPI/ERA-I/';
% % %GPIdir    = '/Users/meganfowler/gp_fuse/TC/GPI_daily/Bootstrap_testResults/';
% GPIdir    = '/Users/meganfowler/gp_fuse/TC/GPI_daily/Bootstrap_testResults/';
% % 
% % %% Test section
% % % GPIdir   = '/Volumes/MyPassport/Data/TCs/GPI/ERA-I/GPI_DailyGPIvars_MJOphase1.nc';  %Old file
% % % %newDir   = '/Users/meganfowler/gp_fuse/TC/GPI_daily/Bootstrap_testResults/GPI_DailyGPIvars_MJOphase_1-OMI_iBoot1.nc'; %Bootstrapped file
% % % newDir   = '~/Desktop/GPI_DailyGPIvars_MJOphase_1-OMI_iBoot1.nc'; %Bootstrapped file
% % % 
% % % 
% % % oldLon = ncread(GPIdir,'lon');
% % % oldLat = ncread(GPIdir,'lat');
% % % ilon = find(oldLon>=60 & oldLon<=200);
% % % ilat = find(oldLat>-5 & oldLat<40);
% % % 
% % % var_old = ncread(GPIdir,'GPI'); %[lon lat time]
% % % var_new = ncread(newDir,'GPI'); %[lon lat time]
% % % var_old(var_old<=-32767) = NaN;
% % % var_new(var_new<=-32767) = NaN;
% % % 
% % % diffVars = squeeze((var_new(:,:,iJune)-var_old(ilon,ilat,iJune)));
% % % diffPct = diffVars./squeeze(var_old(ilon,ilat,iJune));
% % % 
% % % figure;
% % % contourf(lon,lat,isfinite(diffPct)','linecolor','none')
% % % hold on;
% % % plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k');
% % % axis([50 210 0 45]);
% % % colorbar()
% % % title('Diff (%) GPI')
% % % 
% % % 
% % % 
% % % 
% % % % Big differences showing up in potential intensity, which uses SST, SLP,
% % % % and P,T,R profiles. 
% % % newFile = '/Users/meganfowler/gp_fuse/TC/GPI_daily/Bootstrap_DailyGPIvariables/DailyGPIvars_MJOphase_1-OMI_iBoot1.nc';
% % % newFile2 = '/Users/meganfowler/gp_fuse/TC/GPI_daily/Bootstrap_DailyGPIvariables/DailyGPIvars_MJOphase_1-OMI_iBoot2.nc';
% % % oldFile = '/Users/meganfowler/gp_fuse/TC/GPI_daily/DailyGPIvars_MJOphase1.nc'; 
% % % 
% % % var_old = ncread(oldFile,'SST'); %[lon lat time]
% % % var_new = ncread(newFile,'SST'); %[lon lat time]
% % % var_new2 = ncread(newFile2,'SST');
% % % 
% % % diffVars = squeeze((var_new(iJune,:,:,1)-var_old(iJune,ilon,ilat,1)));
% % % diffPct = diffVars./squeeze(var_old(iJune,ilon,ilat,1));
% % % 
% % % figure;
% % % contourf(lon,lat,diffPct','linecolor','none')
% % % hold on;
% % % plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k');
% % % axis([50 210 0 45]);
% % % colorbar()
% % % title('Diff (%) in V (level=1)')
% % 
% % 
% %% Define GPI "climatology" as average of all 8 phases 
% nBoot = 2;
% 
% GPI_mjo   = NaN(nBoot,numel(lon),numel(lat),iNov-iJune+1,8);   %Original GPI for mjo phases 
% termsMJO  = NaN(nBoot,numel(lon),numel(lat),iNov-iJune+1,4,8); 
% 
% for iBoot=1:nBoot
%     tic
%     for iPhase=1:8
%         MJOfile = [GPIdir,'GPI_DailyGPIvars_MJOphase_',num2str(iPhase),'-OMI_iBoot',num2str(iBoot),'.nc'];  %OMI index
%         %MJOfile = [GPIdir,'GPI_DailyGPIvars_MJOphase',num2str(iPhase),'.nc'];   %OMI index
% 
%         absVort_mjo = ncread(MJOfile,'ABSVORT'); %[lon lat time]
%         RH700_mjo   = ncread(MJOfile,'RELHUM');
%         Vmax_mjo    = ncread(MJOfile,'VMAX');
%         Vshear_mjo  = ncread(MJOfile,'VSHEAR');
%         tempGPI_mjo = ncread(MJOfile,'GPI'); 
%         GPI_mjo(iBoot,:,:,:,iPhase) = tempGPI_mjo(:,:,iJune:iNov); 
% 
%         absVort_mjo(absVort_mjo<=-32767) = NaN;
%         RH700_mjo(RH700_mjo<=-32767) = NaN;
%         Vmax_mjo(Vmax_mjo<=-32767) = NaN;
%         Vshear_mjo(Vshear_mjo<=-32767) = NaN;
%         GPI_mjo(GPI_mjo<=-32767) = NaN;
% 
%         term1_mjo = abs((10^5) * absVort_mjo).^(3/2);
%         term2_mjo = (RH700_mjo / 50).^3;
%         term3_mjo = (Vmax_mjo / 70).^3;
%         term4_mjo = (1 + (0.1*Vshear_mjo)).^-2;
%         
%         
%         %Make land mask using term 3 (potential intensity term) 
%         mask = double(isfinite(term3_mjo));
%         mask(mask==0) = NaN;
% 
%         term1_mjo = term1_mjo.*mask;
%         term2_mjo = term2_mjo.*mask;
%         term3_mjo = term3_mjo.*mask;  %Unnecessary but consistent 
%         term4_mjo = term4_mjo.*mask;
%         
%         termsMJO(iBoot,:,:,:,1,iPhase) = term1_mjo(:,:,iJune:iNov);
%         termsMJO(iBoot,:,:,:,2,iPhase) = term2_mjo(:,:,iJune:iNov);
%         termsMJO(iBoot,:,:,:,3,iPhase) = term3_mjo(:,:,iJune:iNov);
%         termsMJO(iBoot,:,:,:,4,iPhase) = term4_mjo(:,:,iJune:iNov);   
%     end
%     fprintf('Done with bootstrap %i \n',iBoot); 
%     toc
% end
% 
% 
% %Take average of all the terms, and then use that average to compute GPI
% termsClim = squeeze(nanmean(termsMJO,6));  %Average over phases 
% 
% % GPI_clim = termsClim(:,:,:,1).*termsClim(:,:,:,2).*termsClim(:,:,:,3).*termsClim(:,:,:,4);
% % GPIratio = (termsMJO(:,:,:,1,1)./termsClim(:,:,:,1)).*...
% %     (termsMJO(:,:,iJune:iNov,2,1)./termsClim(:,:,iJune:iNov,2)).*...
% %     (termsMJO(:,:,iJune:iNov,3,1)./termsClim(:,:,iJune:iNov,3)).*...
% %     (termsMJO(:,:,iJune:iNov,4,1)./termsClim(:,:,iJune:iNov,4));
% % 
% %% Start testing out log-ifying
% 
% %Get difference for each phase from climo average
% climSeasonal = squeeze(nanmean(termsClim(:,:,:,:,:),4));
% 
% for iPhase=1:8
%     iCount=1;
%     for iDay=1:(iNov-iJune+1)
%         diffTerm(:,:,:,iCount,:,iPhase) = squeeze(termsMJO(:,:,:,iDay,:,iPhase))-climSeasonal(:,:,:,:); 
%         iCount = iCount+1;
%     end 
% end
% 
% %As in eq. 5 of Zhao and Li (2018), define coefficients as product of mean
% %   terms not being changed in a given equation 
% alpha1 = termsClim(:,:,:,:,2).*termsClim(:,:,:,:,3).*termsClim(:,:,:,:,4);
% alpha2 = termsClim(:,:,:,:,1).*termsClim(:,:,:,:,3).*termsClim(:,:,:,:,4);
% alpha3 = termsClim(:,:,:,:,1).*termsClim(:,:,:,:,2).*termsClim(:,:,:,:,4);
% alpha4 = termsClim(:,:,:,:,1).*termsClim(:,:,:,:,2).*termsClim(:,:,:,:,3); 
% 
% %Use sum to get at full GPI in each phase of the MJO 
% for iPhase=1:8
%     dGPI_phs(:,:,:,:,iPhase) = (alpha1.*squeeze(diffTerm(:,:,:,:,1,iPhase))) + (alpha2.*squeeze(diffTerm(:,:,:,:,2,iPhase))) + ...
%                              (alpha3.*squeeze(diffTerm(:,:,:,:,3,iPhase))) + (alpha4.*squeeze(diffTerm(:,:,:,:,4,iPhase)));                       
% end
% 
% %Get regional averages of each term 
% lonsLeft  = [100,120,130,140,160];
% lonsRight = [120,130,140,160,180];
% ilats = find(lat>=5 & lat<=20);
% 
% fprintf('Now at line 153 - getting averages for each term/bootstrap/region/phase');
% 
% for iBoot=1:nBoot
%     for iPhase=1:8
%         for iReg=1:5
%            ilons = find(lon>=lonsLeft(iReg) & lon<lonsRight(iReg)); 
% 
%            dGPI_avg(iBoot,iPhase,iReg)     = nanmean(nanmean(nanmean(dGPI_phs(iBoot,ilons,ilats,:,iPhase))));
%            term1_avg(iBoot,iPhase,iReg)    = nanmean(nanmean(nanmean((squeeze(alpha1(iBoot,ilons,ilats,:)).*squeeze(diffTerm(iBoot,ilons,ilats,:,1,iPhase))))));
%            term2_avg(iBoot,iPhase,iReg)    = nanmean(nanmean(nanmean((squeeze(alpha2(iBoot,ilons,ilats,:)).*squeeze(diffTerm(iBoot,ilons,ilats,:,2,iPhase))))));
%            term3_avg(iBoot,iPhase,iReg)    = nanmean(nanmean(nanmean((squeeze(alpha3(iBoot,ilons,ilats,:)).*squeeze(diffTerm(iBoot,ilons,ilats,:,3,iPhase))))));
%            term4_avg(iBoot,iPhase,iReg)    = nanmean(nanmean(nanmean((squeeze(alpha4(iBoot,ilons,ilats,:)).*squeeze(diffTerm(iBoot,ilons,ilats,:,4,iPhase))))));
% 
%            clearvars ilons 
%         end
%     end
% end

%Quick check that the map still looks reasonable for fullGPI phase 1:
%   contourf(lon,lat,nanmean(dGPI_phs(:,:,:,1),3)',-4:0.04:4,'linecolor','none'); 

%% Try plotting 
% 
% iGPI = 1; 
% c = categorical({'100?-120?','120?-130?','130?-140?','140?-160?','160?-180?'});
% 
% phs1 = [term1_avg(iGPI,1),term2_avg(iGPI,1),term3_avg(iGPI,1),term4_avg(iGPI,1),dGPI_avg(iGPI,1);...
%         term1_avg(iGPI,2),term2_avg(iGPI,2),term3_avg(iGPI,2),term4_avg(iGPI,2),dGPI_avg(iGPI,2);...
%         term1_avg(iGPI,3),term2_avg(iGPI,3),term3_avg(iGPI,3),term4_avg(iGPI,3),dGPI_avg(iGPI,3);...
%         term1_avg(iGPI,4),term2_avg(iGPI,4),term3_avg(iGPI,4),term4_avg(iGPI,4),dGPI_avg(iGPI,4);...
%         term1_avg(iGPI,5),term2_avg(iGPI,5),term3_avg(iGPI,5),term4_avg(iGPI,5),dGPI_avg(iGPI,5)];
% 
% figure;
% bar(c,phs1);
% legend({'Vorticity','Humidity','Potential Intensity','Shear','fullGPI'},'location','northwest');
% xlabel('Region','fontsize',20)
% title('Phase 1','fontsize',20);
% ax = gca;
% ax.FontSize = 16;
% 
% % MOVE INTO FOR LOOP %
% figure;
% c = categorical({'100?-120?','120?-130?','130?-140?','140?-160?','160?-180?'});
% 
% for iGPI=1:8    
%     phsRegions = [term1_avg(iGPI,1),term2_avg(iGPI,1),term3_avg(iGPI,1),term4_avg(iGPI,1),dGPI_avg(iGPI,1);...
%             term1_avg(iGPI,2),term2_avg(iGPI,2),term3_avg(iGPI,2),term4_avg(iGPI,2),dGPI_avg(iGPI,2);...
%             term1_avg(iGPI,3),term2_avg(iGPI,3),term3_avg(iGPI,3),term4_avg(iGPI,3),dGPI_avg(iGPI,3);...
%             term1_avg(iGPI,4),term2_avg(iGPI,4),term3_avg(iGPI,4),term4_avg(iGPI,4),dGPI_avg(iGPI,4);...
%             term1_avg(iGPI,5),term2_avg(iGPI,5),term3_avg(iGPI,5),term4_avg(iGPI,5),dGPI_avg(iGPI,5)];
%         
%     figure;
%     %Only add x-axis regions if bottom plot
%     if iGPI<8
%         bar(phsRegions);
%         ax = gca;
%         ax.FontSize = 14;
%         title(['Phase ',num2str(iGPI)],'fontsize',16);
%     elseif iGPI==8
%         bar(c,phsRegions);
%         xlabel('Region','FontSize',16)
%         title(['Phase ',num2str(iGPI)],'fontsize',16);
%         ax = gca;
%         ax.FontSize = 14;
%     end
%     
%     %Only add legend to top plot
%     if iGPI==1
%         legend({'Vorticity','Humidity','Potential Intensity','Shear','fullGPI'},'location','northwest');
%     end
%     
%        
%     clearvars phsRegions
%     %print(['~/Desktop/TestBarChart_Phs',num2str(iGPI),'.pdf'],'-dpdf','-bestfit')
% end

%% Try plotting the mean based on the bootsrap with error bars 

%Took too long to bus between pritchnode and local, so copied over this
%   code to /fast/mdfowler and ran it on the group node - super fast. Saved
%   out termX_avg and dGPI_avg and nBoot. 
load('/Users/meganfowler/gp_fuse/TC/GPI_daily/bootstrap_scripts/BootstrappedGPIterms_100_5to20lat.mat');

term1_avgBootstraps = squeeze(nanmean(term1_avg,1));
term2_avgBootstraps = squeeze(nanmean(term2_avg,1));
term3_avgBootstraps = squeeze(nanmean(term3_avg,1));
term4_avgBootstraps = squeeze(nanmean(term4_avg,1));
dGPI_avgBootstraps  = squeeze(nanmean(dGPI_avg,1)); 

% stdErr_term1 = squeeze(nanstd(term1_avg))/sqrt(nBoot-1);
% stdErr_term2 = squeeze(nanstd(term2_avg))/sqrt(nBoot-1);
% stdErr_term3 = squeeze(nanstd(term3_avg))/sqrt(nBoot-1);
% stdErr_term4 = squeeze(nanstd(term4_avg))/sqrt(nBoot-1);
% stdErr_dGPI  = squeeze(nanstd(dGPI_avg))/sqrt(nBoot-1);

term1_pctile = prctile(term1_avg,[25,75],1);
term2_pctile = prctile(term2_avg,[25,75],1);
term3_pctile = prctile(term3_avg,[25,75],1);
term4_pctile = prctile(term4_avg,[25,75],1);
dGPI_pctile  = prctile(dGPI_avg,[25,75],1);

%% Line plot?
% figure;
% errorbar(1:8,dGPI_avgBootstraps(:,1),dGPI_pctile(1,:,1), dGPI_pctile(2,:,1),'red','LineWidth',.75)
% hold on 
% errorbar(1:8,dGPI_avgBootstraps(:,2),dGPI_pctile(1,:,2), dGPI_pctile(2,:,2),'Color',[0.9290 0.6940 0.1250],'LineWidth',.75)
% errorbar(1:8,dGPI_avgBootstraps(:,3),dGPI_pctile(1,:,3), dGPI_pctile(2,:,3),'Color',[0.4660 0.6740 0.1880],'LineWidth',.75)
% errorbar(1:8,dGPI_avgBootstraps(:,4),dGPI_pctile(1,:,4), dGPI_pctile(2,:,4),'Color',[0 0.4470 0.7410],'LineWidth',.75)
% errorbar(1:8,dGPI_avgBootstraps(:,5),dGPI_pctile(1,:,5), dGPI_pctile(2,:,5),'Color',[0.4940 0.1840 0.5560],'LineWidth',.75)
% legend('100?-120?','120?-130?','130?-140?','140?-160?','160?-180?')
% %Plot thicker main line
% plot(1:8, dGPI_avgBootstraps(:,1), 'red','LineWidth',2,'HandleVisibility','off'); 
% plot(1:8, dGPI_avgBootstraps(:,2), 'Color',[0.9290 0.6940 0.1250],'LineWidth',2,'HandleVisibility','off'); 
% plot(1:8, dGPI_avgBootstraps(:,3), 'Color',[0.4660 0.6740 0.1880],'LineWidth',2,'HandleVisibility','off'); 
% plot(1:8, dGPI_avgBootstraps(:,4), 'Color',[0 0.4470 0.7410],'LineWidth',2,'HandleVisibility','off'); 
% plot(1:8, dGPI_avgBootstraps(:,5), 'Color',[0.4940 0.1840 0.5560],'LineWidth',2,'HandleVisibility','off'); 
% %Plot options
% plot([0.5,8.5],[0 0],'k-','HandleVisibility','off')  %Plot solid black line at 0 
% title('GPI progression','FontSize',20)
% xlabel('MJO Phase','FontSize',16);
% ylabel('GPI anomaly','FontSize',16); 
% xlim([0.5 8.5])

%% Bar chart with just GPI? 
figure; 
regionLabels = {'100?-120?','120?-130?','130?-140?','140?-160?','160?-180?'};

%Plot bars and change colors
hBar = bar(1:8,dGPI_avgBootstraps,'FaceColor','flat');
legend(regionLabels,'fontsize',16)
hBar(1).CData = [1 0 0];
hBar(2).CData = [0.9290 0.6940 0.1250];
hBar(3).CData = [0.4660 0.6740 0.1880];
hBar(4).CData = [0 0.4470 0.7410];
hBar(5).CData = [0.4940 0.1840 0.5560];
%Add error bars
for k1=1:5
    ctr(k1,:) = bsxfun(@plus, hBar(1).XData, [hBar(k1).XOffset]');
    ydt(k1,:) = hBar(k1).YData;
    %hBar(k1).CData = k1; 
end
hold on;
errorbar(ctr',ydt',squeeze(dGPI_pctile(1,:,:)),squeeze(dGPI_pctile(2,:,:)),'Color',[0.4,0.4,0.4],'linestyle','none','LineWidth',.75,'HandleVisibility','off')
%Plot options
title('GPI progression', 'fontsize',24); 
ylabel('GPI anomaly', 'fontsize',20); 
xlabel('MJO Phase (OMI)', 'fontsize',20);
ax=gca;
ax.FontSize=16; 
%Add lines to separate phases
plot([1.5 1.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')
plot([2.5 2.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')
plot([3.5 3.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')
plot([4.5 4.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')
plot([5.5 5.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')
plot([6.5 6.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')
plot([7.5 7.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')


%% Figure with rows as regions instead of phases (more readable?) 

figure;
c            = categorical({'Phase 1', 'Phase 2', 'Phase 3', 'Phase 4', 'Phase 5', 'Phase 6', 'Phase 7', 'Phase 8'});
regionLabels = {'100?-120?','120?-130?','130?-140?','140?-160?','160?-180?'};

for iGPI=1:5   
    phsRegions = [term1_avgBootstraps(1,iGPI),term2_avgBootstraps(1,iGPI),term3_avgBootstraps(1,iGPI),term4_avgBootstraps(1,iGPI),dGPI_avgBootstraps(1,iGPI);...
                  term1_avgBootstraps(2,iGPI),term2_avgBootstraps(2,iGPI),term3_avgBootstraps(2,iGPI),term4_avgBootstraps(2,iGPI),dGPI_avgBootstraps(2,iGPI);...
                  term1_avgBootstraps(3,iGPI),term2_avgBootstraps(3,iGPI),term3_avgBootstraps(3,iGPI),term4_avgBootstraps(3,iGPI),dGPI_avgBootstraps(3,iGPI);...
                  term1_avgBootstraps(4,iGPI),term2_avgBootstraps(4,iGPI),term3_avgBootstraps(4,iGPI),term4_avgBootstraps(4,iGPI),dGPI_avgBootstraps(4,iGPI);...
                  term1_avgBootstraps(5,iGPI),term2_avgBootstraps(5,iGPI),term3_avgBootstraps(5,iGPI),term4_avgBootstraps(5,iGPI),dGPI_avgBootstraps(5,iGPI);...
                  term1_avgBootstraps(6,iGPI),term2_avgBootstraps(6,iGPI),term3_avgBootstraps(6,iGPI),term4_avgBootstraps(6,iGPI),dGPI_avgBootstraps(6,iGPI);...
                  term1_avgBootstraps(7,iGPI),term2_avgBootstraps(7,iGPI),term3_avgBootstraps(7,iGPI),term4_avgBootstraps(7,iGPI),dGPI_avgBootstraps(7,iGPI);...
                  term1_avgBootstraps(8,iGPI),term2_avgBootstraps(8,iGPI),term3_avgBootstraps(8,iGPI),term4_avgBootstraps(8,iGPI),dGPI_avgBootstraps(8,iGPI)];
     

     errHigh   = [term1_pctile(2,1,iGPI),term2_pctile(2,1,iGPI),term3_pctile(2,1,iGPI),term4_pctile(2,1,iGPI),dGPI_pctile(2,1,iGPI);...
                  term1_pctile(2,2,iGPI),term2_pctile(2,2,iGPI),term3_pctile(2,2,iGPI),term4_pctile(2,2,iGPI),dGPI_pctile(2,2,iGPI);...
                  term1_pctile(2,3,iGPI),term2_pctile(2,3,iGPI),term3_pctile(2,3,iGPI),term4_pctile(2,3,iGPI),dGPI_pctile(2,3,iGPI);...
                  term1_pctile(2,4,iGPI),term2_pctile(2,4,iGPI),term3_pctile(2,4,iGPI),term4_pctile(2,4,iGPI),dGPI_pctile(2,4,iGPI);...
                  term1_pctile(2,5,iGPI),term2_pctile(2,5,iGPI),term3_pctile(2,5,iGPI),term4_pctile(2,5,iGPI),dGPI_pctile(2,5,iGPI);...
                  term1_pctile(2,6,iGPI),term2_pctile(2,6,iGPI),term3_pctile(2,6,iGPI),term4_pctile(2,6,iGPI),dGPI_pctile(2,6,iGPI);...
                  term1_pctile(2,7,iGPI),term2_pctile(2,7,iGPI),term3_pctile(2,7,iGPI),term4_pctile(2,7,iGPI),dGPI_pctile(2,7,iGPI);...
                  term1_pctile(2,8,iGPI),term2_pctile(2,8,iGPI),term3_pctile(2,8,iGPI),term4_pctile(2,8,iGPI),dGPI_pctile(2,8,iGPI)];
         
     errLow    = [term1_pctile(1,1,iGPI),term2_pctile(1,1,iGPI),term3_pctile(1,1,iGPI),term4_pctile(1,1,iGPI),dGPI_pctile(1,1,iGPI);...
                  term1_pctile(1,2,iGPI),term2_pctile(1,2,iGPI),term3_pctile(1,2,iGPI),term4_pctile(1,2,iGPI),dGPI_pctile(1,2,iGPI);...
                  term1_pctile(1,3,iGPI),term2_pctile(1,3,iGPI),term3_pctile(1,3,iGPI),term4_pctile(1,3,iGPI),dGPI_pctile(1,3,iGPI);...
                  term1_pctile(1,4,iGPI),term2_pctile(1,4,iGPI),term3_pctile(1,4,iGPI),term4_pctile(1,4,iGPI),dGPI_pctile(1,4,iGPI);...
                  term1_pctile(1,5,iGPI),term2_pctile(1,5,iGPI),term3_pctile(1,5,iGPI),term4_pctile(1,5,iGPI),dGPI_pctile(1,5,iGPI);...
                  term1_pctile(1,6,iGPI),term2_pctile(1,6,iGPI),term3_pctile(1,6,iGPI),term4_pctile(1,6,iGPI),dGPI_pctile(1,6,iGPI);...
                  term1_pctile(1,7,iGPI),term2_pctile(1,7,iGPI),term3_pctile(1,7,iGPI),term4_pctile(1,7,iGPI),dGPI_pctile(1,7,iGPI);...
                  term1_pctile(1,8,iGPI),term2_pctile(1,8,iGPI),term3_pctile(1,8,iGPI),term4_pctile(1,8,iGPI),dGPI_pctile(1,8,iGPI)];
              
%     errRegions = [stdErr_term1(1,iGPI),stdErr_term2(1,iGPI),stdErr_term3(1,iGPI),stdErr_term4(1,iGPI),stdErr_dGPI(1,iGPI);...
%                   stdErr_term1(2,iGPI),stdErr_term2(2,iGPI),stdErr_term3(2,iGPI),stdErr_term4(2,iGPI),stdErr_dGPI(2,iGPI);...
%                   stdErr_term1(3,iGPI),stdErr_term2(3,iGPI),stdErr_term3(3,iGPI),stdErr_term4(3,iGPI),stdErr_dGPI(3,iGPI);...
%                   stdErr_term1(4,iGPI),stdErr_term2(4,iGPI),stdErr_term3(4,iGPI),stdErr_term4(4,iGPI),stdErr_dGPI(4,iGPI);...
%                   stdErr_term1(5,iGPI),stdErr_term2(5,iGPI),stdErr_term3(5,iGPI),stdErr_term4(5,iGPI),stdErr_dGPI(5,iGPI);...
%                   stdErr_term1(6,iGPI),stdErr_term2(6,iGPI),stdErr_term3(6,iGPI),stdErr_term4(6,iGPI),stdErr_dGPI(6,iGPI);...
%                   stdErr_term1(7,iGPI),stdErr_term2(7,iGPI),stdErr_term3(7,iGPI),stdErr_term4(7,iGPI),stdErr_dGPI(7,iGPI);...
%                   stdErr_term1(8,iGPI),stdErr_term2(8,iGPI),stdErr_term3(8,iGPI),stdErr_term4(8,iGPI),stdErr_dGPI(8,iGPI)];

        
    %Only add x-axis regions if bottom plot
    subplot(3,2,iGPI)     
    ctrs = 1:8;
    
    hBar = bar(ctrs,phsRegions,1);
    for k1=1:size(phsRegions,2)
        ctr(k1,:) = bsxfun(@plus, hBar(1).XData, [hBar(k1).XOffset]');
        ydt(k1,:) = hBar(k1).YData;
    end
    hold on
    %errorbar(ctr',ydt',2.*errRegions,'k','linestyle','none','LineWidth',1)
    errorbar(ctr',ydt',errLow,errHigh,'Color',[0.4,0.4,0.4],'linestyle','none','LineWidth',.75)

    ylim([-3, 3]);
    xlabel('MJO Phase [OMI]','FontSize',16)
    title(['Region: ',regionLabels{iGPI}],'fontsize',16);
    ax = gca;
    ax.FontSize = 14;
    %Add lines to separate phases
    plot([1.5 1.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')
    plot([2.5 2.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')
    plot([3.5 3.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')
    plot([4.5 4.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')
    plot([5.5 5.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')
    plot([6.5 6.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')
    plot([7.5 7.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')
   

    
    %Only add legend to top plot
    if iGPI==5
        legend({'Vorticity','Humidity','Potential Intensity','Shear','fullGPI'},'location','southeast','fontsize',10);
    end
    
       
    clearvars phsRegions
    %print(['~/Desktop/TestBarChart_Phs',num2str(iGPI),'.pdf'],'-dpdf','-bestfit')
end

fig = gcf;
fig.PaperPositionMode = 'auto';
%print('~/Desktop/TestBarChart_withInterquartileErrBar','-depsc')
print('/Users/meganfowler/Documents/Irvine/TCs/Figures/GPI-BarChart_noMap_5to20lat&landMask','-depsc');
% 
% % % fig.PaperUnits = 'inches';
% % % fig.PaperPosition=[0 0 8 12];
% % % print('5by3DimensionsFigure','-dpng','-r0')
