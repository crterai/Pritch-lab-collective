% Create bar chart for GPI term contributions to full GPI, as plotted in
% map form on the AMS poster and most of my talks. 
%
% 2019-03-05 

%% Read in lat/lon from extra file
SSTfile  = '/Volumes/MyPassport/Data/TCs/globalSST_1983-2013_dailyAvg.nc';
lon = ncread(SSTfile,'longitude');
lat = ncread(SSTfile,'latitude');

coast = load('coast.mat');

%% Read in GPI variables and carry out decomposition 
GPIdir   = '/Volumes/MyPassport/Data/TCs/GPI/ERA-I/';

%% Define GPI "climatology" as average of all 8 phases 

GPI_mjo   = NaN(numel(lon),numel(lat),365,8);   %Original GPI for mjo phases 
termsMJO  = NaN(numel(lon),numel(lat),365,4,8); 

for iPhase=1:8
    MJOfile = [GPIdir,'GPI_DailyGPIvars_MJOphase_',num2str(iPhase),'-RMM.nc'];  %RMM index
    %MJOfile = [GPIdir,'GPI_DailyGPIvars_MJOphase',num2str(iPhase),'.nc'];   %OMI index
    
    absVort_mjo = ncread(MJOfile,'ABSVORT'); %[lon lat time]
    RH700_mjo   = ncread(MJOfile,'RELHUM');
    Vmax_mjo    = ncread(MJOfile,'VMAX');
    Vshear_mjo  = ncread(MJOfile,'VSHEAR');
    GPI_mjo(:,:,:,iPhase) = ncread(MJOfile,'GPI'); 

    absVort_mjo(absVort_mjo<=-32767) = NaN;
    RH700_mjo(RH700_mjo<=-32767) = NaN;
    Vmax_mjo(Vmax_mjo<=-32767) = NaN;
    Vshear_mjo(Vshear_mjo<=-32767) = NaN;
    GPI_mjo(GPI_mjo<=-32767) = NaN;
       
    term1_mjo = abs((10^5) * absVort_mjo).^(3/2);
    term2_mjo = (RH700_mjo / 50).^3;
    term3_mjo = (Vmax_mjo / 70).^3;
    term4_mjo = (1 + (0.1*Vshear_mjo)).^-2; 

    termsMJO(:,:,:,1,iPhase) = term1_mjo;
    termsMJO(:,:,:,2,iPhase) = term2_mjo;
    termsMJO(:,:,:,3,iPhase) = term3_mjo;
    termsMJO(:,:,:,4,iPhase) = term4_mjo;   
end

%Take average of all the terms, and then use that average to compute GPI
termsClim = squeeze(nanmean(termsMJO,5));
iJune = 152; %Day of year that June starts
iNov  = 334; %Day of year that November ends 


GPI_clim = termsClim(:,:,:,1).*termsClim(:,:,:,2).*termsClim(:,:,:,3).*termsClim(:,:,:,4);
GPIratio = (termsMJO(:,:,iJune:iNov,1,1)./termsClim(:,:,iJune:iNov,1)).*...
    (termsMJO(:,:,iJune:iNov,2,1)./termsClim(:,:,iJune:iNov,2)).*...
    (termsMJO(:,:,iJune:iNov,3,1)./termsClim(:,:,iJune:iNov,3)).*...
    (termsMJO(:,:,iJune:iNov,4,1)./termsClim(:,:,iJune:iNov,4));

%% Start testing out log-ifying

%Get difference for each phase from climo average
climSeasonal = squeeze(nanmean(termsClim(:,:,iJune:iNov,:),3));

for iPhase=1:8
    iCount=1;
    for iDay=iJune:iNov
        diffTerm(:,:,iCount,:,iPhase) = squeeze(termsMJO(:,:,iDay,:,iPhase))-climSeasonal(:,:,:); 
        iCount = iCount+1;
    end 
end

%As in eq. 5 of Zhao and Li (2018), define coefficients as product of mean
%   terms not being changed in a given equation 
alpha1 = termsClim(:,:,iJune:iNov,2).*termsClim(:,:,iJune:iNov,3).*termsClim(:,:,iJune:iNov,4);
alpha2 = termsClim(:,:,iJune:iNov,1).*termsClim(:,:,iJune:iNov,3).*termsClim(:,:,iJune:iNov,4);
alpha3 = termsClim(:,:,iJune:iNov,1).*termsClim(:,:,iJune:iNov,2).*termsClim(:,:,iJune:iNov,4);
alpha4 = termsClim(:,:,iJune:iNov,1).*termsClim(:,:,iJune:iNov,2).*termsClim(:,:,iJune:iNov,3); 

%Use sum to get at full GPI in each phase of the MJO 
for iPhase=1:8
    dGPI_phs(:,:,:,iPhase) = (alpha1.*squeeze(diffTerm(:,:,:,1,iPhase))) + (alpha2.*squeeze(diffTerm(:,:,:,2,iPhase))) + ...
                             (alpha3.*squeeze(diffTerm(:,:,:,3,iPhase))) + (alpha4.*squeeze(diffTerm(:,:,:,4,iPhase)));                       
end

%Get regional averages of each term 
lonsLeft  = [100,120,130,140,160];
lonsRight = [120,130,140,160,180];
ilats = find(lat>=0 & lat<=30);

for iPhase=1:8
    for iReg=1:5
       ilons = find(lon>=lonsLeft(iReg) & lon<lonsRight(iReg)); 

       dGPI_avg(iPhase,iReg)     = nanmean(nanmean(nanmean(dGPI_phs(ilons,ilats,:,iPhase))));
       term1_avg(iPhase,iReg)    = nanmean(nanmean(nanmean((alpha1(ilons,ilats,:).*squeeze(diffTerm(ilons,ilats,:,1,iPhase))))));
       term2_avg(iPhase,iReg)    = nanmean(nanmean(nanmean((alpha2(ilons,ilats,:).*squeeze(diffTerm(ilons,ilats,:,2,iPhase))))));
       term3_avg(iPhase,iReg)    = nanmean(nanmean(nanmean((alpha3(ilons,ilats,:).*squeeze(diffTerm(ilons,ilats,:,3,iPhase))))));
       term4_avg(iPhase,iReg)    = nanmean(nanmean(nanmean((alpha4(ilons,ilats,:).*squeeze(diffTerm(ilons,ilats,:,4,iPhase))))));

       clearvars ilons 
    end
end

%Quick check that the map still looks reasonable for fullGPI phase 1:
%   contourf(lon,lat,nanmean(dGPI_phs(:,:,:,1),3)',-4:0.04:4,'linecolor','none'); 

%% Bar chart with just GPI? 
figure; 
regionLabels = {'100?-120?','120?-130?','130?-140?','140?-160?','160?-180?'};

%Plot bars and change colors
hBar = bar(1:8,dGPI_avg,'FaceColor','flat');
legend(regionLabels,'fontsize',10)
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
%Plot options
title('GPI progression', 'fontsize',20); 
ylabel('GPI anomaly', 'fontsize',16); 
xlabel('MJO Phase', 'fontsize',16);
%Add lines to separate phases
plot([1.5 1.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')
plot([2.5 2.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')
plot([3.5 3.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')
plot([4.5 4.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')
plot([5.5 5.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')
plot([6.5 6.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')
plot([7.5 7.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')

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

%% Figure with rows as regions instead of phases (more readable?) 

figure;
c            = categorical({'Phase 1', 'Phase 2', 'Phase 3', 'Phase 4', 'Phase 5', 'Phase 6', 'Phase 7', 'Phase 8'});
regionLabels = {'100?-120?','120?-130?','130?-140?','140?-160?','160?-180?'};

for iGPI=1:5   
    phsRegions = [term1_avg(1,iGPI),term2_avg(1,iGPI),term3_avg(1,iGPI),term4_avg(1,iGPI),dGPI_avg(1,iGPI);...
                  term1_avg(2,iGPI),term2_avg(2,iGPI),term3_avg(2,iGPI),term4_avg(2,iGPI),dGPI_avg(2,iGPI);...
                  term1_avg(3,iGPI),term2_avg(3,iGPI),term3_avg(3,iGPI),term4_avg(3,iGPI),dGPI_avg(3,iGPI);...
                  term1_avg(4,iGPI),term2_avg(4,iGPI),term3_avg(4,iGPI),term4_avg(4,iGPI),dGPI_avg(4,iGPI);...
                  term1_avg(5,iGPI),term2_avg(5,iGPI),term3_avg(5,iGPI),term4_avg(5,iGPI),dGPI_avg(5,iGPI);...
                  term1_avg(6,iGPI),term2_avg(6,iGPI),term3_avg(6,iGPI),term4_avg(6,iGPI),dGPI_avg(6,iGPI);...
                  term1_avg(7,iGPI),term2_avg(7,iGPI),term3_avg(7,iGPI),term4_avg(7,iGPI),dGPI_avg(7,iGPI);...
                  term1_avg(8,iGPI),term2_avg(8,iGPI),term3_avg(8,iGPI),term4_avg(8,iGPI),dGPI_avg(8,iGPI)];
        
    %Only add x-axis regions if bottom plot
    subplot(3,2,iGPI)     
    bar(phsRegions);
    ylim([-1.5, 1.5]);
    xlabel('MJO Phase [RMM]','FontSize',16)
    title(['Region: ',regionLabels{iGPI}],'fontsize',16);
    ax = gca;
    ax.FontSize = 14;
    
    
    %Only add legend to top plot
    if iGPI==5
        legend({'Vorticity','Humidity','Potential Intensity','Shear','fullGPI'},'location','southeast','fontsize',12);
    end
    
       
    clearvars phsRegions
    %print(['~/Desktop/TestBarChart_Phs',num2str(iGPI),'.pdf'],'-dpdf','-bestfit')
end

%fig = gcf;
%fig.PaperPositionMode = 'auto';
%print('~/Desktop/TestBarChart_ByRegion_RMM','-depsc')

% % % fig.PaperUnits = 'inches';
% % % fig.PaperPosition=[0 0 8 12];
% % % print('5by3DimensionsFigure','-dpng','-r0')
