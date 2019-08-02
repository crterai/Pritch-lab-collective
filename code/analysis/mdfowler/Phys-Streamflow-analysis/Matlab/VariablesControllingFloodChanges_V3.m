% PROJECT: Flooding physiology
%
% Determine changes in other variables from the climate change simulations
% that may be acting to cause changes in flood return period. Particularly,
% focus on regions where the sign of the flood change is in disagreement
% with the sign of precipitation changes. 
%
% Megan D. Fowler, 2017-11-16
% Updated on 2018-02-02 to base conclusions on FloodDriverMask_CESMres,
% which separates out which regions are driven by both, phys, or just rad
% changes. Also attempts to summarize information as graph instead of map.
%
% Updated on 2019-01-23 to add 99th percentile rain rate and sfc
% temperature to the tables in the Supplement. 
% 

coast = load('coast.mat');

%Full case should be second (but moved to 4 for consistency with previous) 
caseNames = {'/Users/meganfowler/gp_fuse/Flooding-physiology/RelatedCESMvariables/cesm1_0_6.1850_prei.1deg.001.clm2.MonthlyConcat_Meg.nc'...
             '/Users/meganfowler/gp_fuse/Flooding-physiology/RelatedCESMvariables/cesm1_0_6.1850_4xco2_fdbgb.1deg.002.clm2.MonthlyConcat_Meg.nc'...
             '/Users/meganfowler/gp_fuse/Flooding-physiology/RelatedCESMvariables/cesm1_0_6.1850_4xco2_fixgb.1deg.001.clm2.MonthlyConcat_Meg.nc'...
             '/Users/meganfowler/gp_fuse/Flooding-physiology/RelatedCESMvariables/cesm1_0_6.1850_4xco2_fulgb.1deg.002.clm2.MonthlyConcat_Meg.nc'};

dailyRainFiles = {'/Users/meganfowler/gp_fuse/Flooding-physiology/RelatedCESMvariables/cesm1_0_6.1850_prei.1deg.001.cam2.h1.DailyRainfall.nc'...
                  '/Users/meganfowler/gp_fuse/Flooding-physiology/RelatedCESMvariables/cesm1_0_6.1850_4xco2_fdbgb.1deg.002.cam2.h1.DailyRainfall.nc'...
                  '/Users/meganfowler/gp_fuse/Flooding-physiology/RelatedCESMvariables/cesm1_0_6.1850_4xco2_fixgb.1deg.001.cam2.h1.DailyRainfall.nc'...
                  '/Users/meganfowler/gp_fuse/Flooding-physiology/RelatedCESMvariables/cesm1_0_6.1850_4xco2_fulgb.1deg.002.cam2.h1.DailyRainfall.nc'};          
         
              
%Read in dimension sizes for files 
lon = ncread(caseNames{1},'lon');
lat = ncread(caseNames{1},'lat');

% Read in all variables - VERY time consuming
for iCase=1:4
    tic
    
    sfcT(iCase,:,:,:) = ncread(caseNames{iCase},'TSA');    %2m Sfc Temp [K]
    prect(iCase,:,:,:) = ncread(dailyRainFiles{iCase},'PRECT'); %[m/s]
    timeT(iCase,:)  = nctimenoleap(caseNames{iCase});
    timeP(iCase,:)  = nctimenoleap(dailyRainFiles{iCase}); 
  
    toc 
    fprintf('Finished with case %d \n', iCase)
end

prect = prect.*86400.*1000; %Convert prect from m/s to mm/day

%% If want to clear everything but original variables (i.e., changing case or driver region) 
clearvars -except caseNames dailyRainFiles lon lat sfcT prect timeT timeP


%%  Compute degrees of freedom 
vecTimeT = datevec(double(timeT(1,:)));   %CTRL and PHYS start in same month (Jan)
monthsT  = vecTimeT(:,2); 

% !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
%       EDIT HERE TO CHANGE CASE
caseChoice = 4; %Set to 2 for RAD // and 3 for PHYS // and 4 for FULL
% !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

%Reshape into annual record 
iYearStarts = 1:365:10950;
for i=1:length(iYearStarts)
   
   if i~=30
        prectYearly(:,:,:,i,:) = prect(:,:,:,iYearStarts(i):iYearStarts(i+1)-1);
        timeYearly(i,:) = timeP(1,iYearStarts(i):iYearStarts(i+1)-1);
        %timeVec  = datevec(double(timeYearly(i,:)));
        %monthsYearly(i,1:365)   = timeVec(:,2); 
   elseif i==30
       prectYearly(:,:,:,i,:)  = prect(:,:,:,iYearStarts(i):end);
       timeYearly(i,:) = timeP(1,iYearStarts(i):end);
       %timeVec  = datevec(double(timeYearly(i,:)));
       %monthsYearly(i,1:365)   = timeVec(:,2); 
   end
   
end
vecTimeP = datevec(double(timeYearly(1,:)));   %CTRL and PHYS start in same month (Jan)
monthsP  = vecTimeP(:,2); 


% !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
%   EDIT HERE & IN LOOP TO CHOOSE DIFFERENTLY DRIVEN REGIONS 
% !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
% Compute reduction for DOF within each region -- Moderately time consuming

% Lon/Lat for BOTH driven regions - [SE US, SE Asia, Central Australia] 
% lonStart = [264;  95; 136];
% lonEnd   = [276; 109; 146];
% latStart = [ 29;  12; -31];
% latEnd   = [ 38;  25; -18];

% Lon/Lat for PHYS driven regions - [West Amazon, Sahel, SE Australia]
lonStart = [279;  0; 138];    %Handle Sahel (iReg=2) specially since it has two regions
lonEnd   = [302; 33; 154];
latStart = [-17;  4; -42];
latEnd   = [4;   15; -26];

% Lon/Lat for RAD driven regions - [India, Horn of Africa]
% lonStart = [55; 12];
% lonEnd   = [100; 52]; 
% latStart = [18; -6];
% latEnd   = [36; 12]; 

for ireg = 1:numel(lonStart)
    %% Case for RAD regions
%     ilon = find(lon>=(lonStart(ireg)) & lon<=(lonEnd(ireg)));
%     ilat = find(lat>=(latStart(ireg)) & lat<=(latEnd(ireg)));        
%     
%     if ireg==1
%         imonthT = find(monthsT>=7 & monthsT<=11); %Jun-Oct indexed as 7-11
%         imonthP = find(monthsP>=6 & monthsP<=10);  %Jun-Oct is indexed as 6-10 (daily data)
%     elseif ireg==2
%         imonthT = find(monthsT>=9 | monthsT==1); %Aug-Dec indexed as 9-12 and 1
%         imonthP = find(monthsP>=8 & monthsP<=12);  %Agu-Dec is indexed as 3-8 (daily data)
%    end
     
    %% Case for PHYS regions  
    if ireg==2 
        ilon = find(lon>=350 & lon<=360 | lon>=0 & lon<=33);
        ilat = find(lat>=(latStart(ireg)) & lat<=(latEnd(ireg)));
    else
        ilon = find(lon>=(lonStart(ireg)) & lon<=(lonEnd(ireg)));
        ilat = find(lat>=(latStart(ireg)) & lat<=(latEnd(ireg)));        
    end
    
    if ireg==1
        imonthT = find(monthsT>=1 & monthsT<=7); %Dec-Jun is indexed as 1-7
        imonthP = find(monthsP>=12 | monthsP<=6);  %Dec-Jun indexed as 12 or 1-6 (daily)
    elseif ireg==2
        imonthT = find(monthsT>=6 & monthsT<=12); %May-Nov is indexed as 6-12
        imonthP = find(monthsP>=5 & monthsP<=11); %May-Nov is indexed as 5-11 (daily)
    elseif ireg==3
        imonthT = find(monthsT==12 | monthsT<=6); %Nov-May is indexed 12 and 1-6
        imonthP = find(monthsP>=11 | monthsP<=5); %Nov-May is indexed as 11-12 and 1-5 (daily) 
    end

    %% Case for BOTH regions
%     ilon = find(lon>=(lonStart(ireg)) & lon<=(lonEnd(ireg)));
%     ilat = find(lat>=(latStart(ireg)) & lat<=(latEnd(ireg)));
% 
%     if ireg==1
%         imonthT = find(monthsT>=4 & monthsT<=9);  %Mar-Aug is indexed as 4-9
%         imonthP = find(monthsP>=3 & monthsP<=8);  %Mar-Aug is indexed as 3-8 (daily data)
%     elseif ireg==2
%         imonthT = find(monthsT>=7 & monthsT<=11); %Jun-Oct is indexed as 7-11
%         imonthP = find(monthsP>=6 & monthsP<=10); %Jun-Oct is indexed as 6-10 (daily data)
%     elseif ireg==3
%         imonthT = find(monthsT>=1 & monthsT<=5); %Dec-Apr is indexed as 1-5
%         imonthP = find(monthsP>=12 | monthsP<=4); %Dec-Apr is indexed as 12 and 1-4
%     end
%%    
    % -- Isolate months and take differences -- %
    
    %99th percentile precip 
    rainReg   = squeeze(prectYearly(caseChoice,ilon,ilat,:,imonthP)); 
    rainCtrl  = squeeze(prectYearly(1,ilon,ilat,:,imonthP));
    %Set infinity values to NaNs 
    rainCtrl(~isfinite(rainCtrl)) = NaN;
    rainReg(~isfinite(rainReg)) = NaN; 
    %Take percentile 
    rain99_case  = prctile(rainReg,99,4);
    rain99_ctrl = prctile(rainCtrl,99,4); 
    
    diff99   = rain99_case - rain99_ctrl; 
    
    %Surface temperature 
    diffTemp =  squeeze(sfcT(caseChoice,ilon,ilat,imonthT) - sfcT(caseChoice,ilon,ilat,imonthT));
           
    %Compute effective ratio for DOF 
    [ratio_eff_pr99, Ratio_per_dim_pr99] = Eff_spatial_DOF_ratio(diff99);
    [ratio_eff_sfcT, Ratio_per_dim_sfcT] = Eff_spatial_DOF_ratio(diffTemp); 
    
    ratio_eff_full(ireg,:) = [ratio_eff_pr99; ratio_eff_sfcT];
        
    clearvars ilon ilat imonthP imonthT diff99 diffTemp rainReg rainCtrl rain99_case rain99_ctrl
    fprintf('Finished with region %d of %d \n', ireg, numel(lonStart));
end

%% 
%Load in the mask
load('/Users/meganfowler/gp_fuse/Flooding-physiology/MatlabData/pyMask_CESMres_95sig.mat'); %1=Both; 2=Phys; 3=Rad; 
iBoth = find(physMaskNew==1); 
iPhys = find(physMaskNew==2);
iRad  = find(physMaskNew==3); 

%Create a binary mask for regions driven by both rad and phys 
maskBoth = NaN(size(physMaskNew));

% !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
%    EDIT HERE (WITHIN THE PARANTHESES) TO CHANGE DRIVERS 
maskBoth(iPhys) = 1; % <- This line needs to be changed according to different runs (iBoth for multiply stressed regions,iPhys for phys regions, etc.)
% !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

%Expand mask to 30 time samples
for t=1:30
    maskBoth_expand30(:,:,t) = maskBoth(:,:);    
end

areaCESM = ncread(caseNames{1},'area');
areaCESM = areaCESM.*maskBoth; %Mask out areas not being considered 

% !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
%   EDIT HERE & IN LOOP TO CHOOSE DIFFERENTLY DRIVEN REGIONS 
% !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

% ---- Create regional bar charts --- 
% Define regions and refine data to that location/season for BOTH forced
% Regions go in the following order:
%   US, SE Asia, Central Australia
% lonStart = [264;  95; 136];
% lonEnd   = [276; 109; 146];
% latStart = [ 29;  12; -31];
% latEnd   = [ 38;  25; -18];
 
% % Lon/Lat for PHYS driven regions - [West Amazon, Sahel, SE Australia]
lonStart = [279;  0; 138];    %Handle Sahel (iReg=2) specially since it has two regions
lonEnd   = [302; 33; 154];
latStart = [-17;  4; -42];
latEnd   = [4;   15; -26];

% Lon/Lat for RAD driven regions - [India, Horn of Africa]
% lonStart = [55; 12];
% lonEnd   = [100; 52]; 
% latStart = [18; -6];
% latEnd   = [36; 12]; 

for ireg = 1:numel(lonStart)
    %% Define restriction on space/time: RAD REGIONS
%     ilon = find(lon>=(lonStart(ireg)) & lon<=(lonEnd(ireg)));
%     ilat = find(lat>=(latStart(ireg)) & lat<=(latEnd(ireg)));        
%     areaReg = areaCESM(ilon,ilat);
%     
%     if ireg==1
%         imonthT = find(monthsT>=7 & monthsT<=11); %Jun-Oct indexed as 7-11
%         imonthP = find(monthsP>=6 & monthsP<=10);  %Jun-Oct is indexed as 6-10 (daily data)
%     elseif ireg==2
%         imonthT = find(monthsT>=9 | monthsT==1); %Aug-Dec indexed as 9-12 and 1
%         imonthP = find(monthsP>=8 & monthsP<=12);  %Agu-Dec is indexed as 3-8 (daily data)
%    end
        
     %% Define restriction on space/time : PHYS REGIONS
    if ireg==2 
        ilon = find(lon>=350 & lon<=360 | lon>=0 & lon<=33);
        ilat = find(lat>=(latStart(ireg)) & lat<=(latEnd(ireg)));
    else
        ilon = find(lon>=(lonStart(ireg)) & lon<=(lonEnd(ireg)));
        ilat = find(lat>=(latStart(ireg)) & lat<=(latEnd(ireg)));        
    end
    areaReg = areaCESM(ilon,ilat);
    
    if ireg==1
        imonthT = find(monthsT>=1 & monthsT<=7);   %Dec-Jun is indexed as 1-7
        imonthP = find(monthsP>=12 | monthsP<=6);  %Dec-Jun indexed as 12 or 1-6 (daily)
    elseif ireg==2
        imonthT = find(monthsT>=6 & monthsT<=12); %May-Nov is indexed as 6-12
        imonthP = find(monthsP>=5 & monthsP<=11); %May-Nov is indexed as 5-11 (daily)
    elseif ireg==3
        imonthT = find(monthsT==12 | monthsT<=6); %Nov-May is indexed 12 and 1-6
        imonthP = find(monthsP>=11 | monthsP<=5); %Nov-May is indexed as 11-12 and 1-5 (daily) 
    end
   
    %% Define restriction on space/time : BOTH REGIONS
%     ilon = find(lon>=(lonStart(ireg)) & lon<=(lonEnd(ireg)));
%     ilat = find(lat>=(latStart(ireg)) & lat<=(latEnd(ireg)));
%     areaReg = areaCESM(ilon,ilat);
%     
%     if ireg==1
%         imonthT = find(monthsT>=4 & monthsT<=9);  %Mar-Aug is indexed as 4-9
%         imonthP = find(monthsP>=3 & monthsP<=8);  %Mar-Aug is indexed as 3-8 (daily data)
%     elseif ireg==2
%         imonthT = find(monthsT>=7 & monthsT<=11); %Jun-Oct is indexed as 7-11
%         imonthP = find(monthsP>=6 & monthsP<=10); %Jun-Oct is indexed as 6-10 (daily data)
%     elseif ireg==3
%         imonthT = find(monthsT>=1 & monthsT<=5); %Dec-Apr is indexed as 1-5
%         imonthP = find(monthsP>=12 | monthsP<=4); %Dec-Apr is indexed as 12 and 1-4
%     end

%  
    %% Regional raw variables, averaged over time, with mask applied
     
    %99th percentile precip 
    rainPhys   = squeeze(prectYearly(3,ilon,ilat,:,imonthP)); %Isolate case, region, and season 
    rainRad    = squeeze(prectYearly(2,ilon,ilat,:,imonthP)); 
    rainFull   = squeeze(prectYearly(4,ilon,ilat,:,imonthP)); 
    rainCtrl   = squeeze(prectYearly(1,ilon,ilat,:,imonthP));
    
    %Set infinity values to NaNs 
    rainCtrl(~isfinite(rainCtrl)) = NaN;
    rainPhys(~isfinite(rainPhys)) = NaN; 
    rainRad(~isfinite(rainRad))   = NaN; 
    rainFull(~isfinite(rainFull)) = NaN; 
    
    rain99phys  = prctile(rainPhys,99,4); %Take 99th percentile precip rate for each year 
    rain99rad   = prctile(rainRad,99,4); 
    rain99full  = prctile(rainFull,99,4);
    rain99ctrl  = prctile(rainCtrl,99,4); 
    
    diffRain_phys = (rain99phys - rain99ctrl).*maskBoth_expand30(ilon,ilat,:); %Define regional difference (with mask applied) 
    diffRain_rad  = (rain99rad  - rain99ctrl).*maskBoth_expand30(ilon,ilat,:);
    diffRain_full = (rain99full - rain99ctrl).*maskBoth_expand30(ilon,ilat,:); 
    
    tempPhys    = squeeze(sfcT(3,ilon,ilat,imonthT));  %Isolate case, region, and season 
    tempRad     = squeeze(sfcT(2,ilon,ilat,imonthT)); 
    tempFull    = squeeze(sfcT(4,ilon,ilat,imonthT)); 
    tempCtrl    = squeeze(sfcT(1,ilon,ilat,imonthT)); 
    
    %Expand mask to nMonths time samples for sfc temperature 
    for t=1:length(imonthT)
        maskBoth_expand(:,:,t) = maskBoth(:,:);    
    end
    diffTemp_phys = (tempPhys - tempCtrl).*maskBoth_expand(ilon,ilat,:); %Define regional difference (with mask applied) 
    diffTemp_rad  = (tempRad  - tempCtrl).*maskBoth_expand(ilon,ilat,:);
    diffTemp_full = (tempFull - tempCtrl).*maskBoth_expand(ilon,ilat,:);
      
    %Area weighting 
    % -- 99th percentile precip 
    rain99phys_avg = squeeze(nanmean(rain99phys,3)).*maskBoth(ilon,ilat);
    rain99rad_avg  = squeeze(nanmean(rain99rad,3)).*maskBoth(ilon,ilat);
    rain99full_avg = squeeze(nanmean(rain99full,3)).*maskBoth(ilon,ilat);
    rain99ctrl_avg = squeeze(nanmean(rain99ctrl,3)).*maskBoth(ilon,ilat);
    
    wgtRainCtrl = rain99ctrl_avg.*areaReg; 
    wgtRainPhys = rain99phys_avg.*areaReg;
    wgtRainRad  = rain99rad_avg.*areaReg; 
    wgtRainFull = rain99full_avg.*areaReg;
    rainCtrlRegWgtAvg   = squeeze(nansum(nansum(wgtRainCtrl)))./squeeze(nansum(nansum(areaReg)));
    rainPhysRegWgtAvg   = squeeze(nansum(nansum(wgtRainPhys)))./squeeze(nansum(nansum(areaReg)));
    rainRadRegWgtAvg    = squeeze(nansum(nansum(wgtRainRad)))./squeeze(nansum(nansum(areaReg)));
    rainFullRegWgtAvg   = squeeze(nansum(nansum(wgtRainFull)))./squeeze(nansum(nansum(areaReg)));
    
    % -- Surface temperature
    tempPhys_avg = squeeze(nanmean(tempPhys,3)).*maskBoth(ilon,ilat); 
    tempRad_avg  = squeeze(nanmean(tempRad,3)).*maskBoth(ilon,ilat);
    tempFull_avg = squeeze(nanmean(tempFull,3)).*maskBoth(ilon,ilat);
    tempCtrl_avg = squeeze(nanmean(tempCtrl,3)).*maskBoth(ilon,ilat);
    
    wgtTempCtrl = tempCtrl_avg.*areaReg;
    wgtTempPhys = tempPhys_avg.*areaReg;
    wgtTempRad  = tempRad_avg.*areaReg;
    wgtTempFull = tempFull_avg.*areaReg;
    tempCtrlRegWgtAvg = squeeze(nansum(nansum(wgtTempCtrl)))./squeeze(nansum(nansum(areaReg)));
    tempPhysRegWgtAvg = squeeze(nansum(nansum(wgtTempPhys)))./squeeze(nansum(nansum(areaReg)));
    tempRadRegWgtAvg  = squeeze(nansum(nansum(wgtTempRad)))./squeeze(nansum(nansum(areaReg)));
    tempFullRegWgtAvg = squeeze(nansum(nansum(wgtTempFull)))./squeeze(nansum(nansum(areaReg)));
        
    %Save all variables
    regionalCtrl(ireg,:) = [rainCtrlRegWgtAvg, tempCtrlRegWgtAvg];
    regionalPhys(ireg,:) = [rainPhysRegWgtAvg, tempPhysRegWgtAvg];
    regionalRad(ireg,:)  = [rainRadRegWgtAvg, tempRadRegWgtAvg];
    regionalFull(ireg,:) = [rainFullRegWgtAvg, tempFullRegWgtAvg];
    
                        
    %% Regional differences

    %Time average 
    rainRegPhys   = squeeze(nanmean(diffRain_phys,3));
    rainRegRad    = squeeze(nanmean(diffRain_rad,3)); 
    rainRegFull   = squeeze(nanmean(diffRain_full,3)); 
    
    tempRegPhys   = squeeze(nanmean(diffTemp_phys,3));
    tempRegRad    = squeeze(nanmean(diffTemp_rad,3));
    tempRegFull   = squeeze(nanmean(diffTemp_full,3)); 
    
    % -- Area weighting -- %
    wgtRainReg_phys = rainRegPhys.*areaReg; 
    wgtRainReg_rad  = rainRegRad.*areaReg; 
    wgtRainReg_full = rainRegFull.*areaReg; 
    
    wgtTempReg_phys = tempRegPhys.*areaReg;
    wgtTempReg_rad  = tempRegRad.*areaReg;
    wgtTempReg_full = tempRegFull.*areaReg; 
 
    %Sum together and divide by sum of area to get area-weighted average
    rainRegWgtAvg_phys   = squeeze(nansum(nansum(wgtRainReg_phys)))./squeeze(nansum(nansum(areaReg)));
    rainRegWgtAvg_rad    = squeeze(nansum(nansum(wgtRainReg_rad)))./squeeze(nansum(nansum(areaReg)));
    rainRegWgtAvg_full   = squeeze(nansum(nansum(wgtRainReg_full)))./squeeze(nansum(nansum(areaReg)));
    
    tempRegWgtAvg_phys   = squeeze(nansum(nansum(wgtTempReg_phys)))./squeeze(nansum(nansum(areaReg)));
    tempRegWgtAvg_rad    = squeeze(nansum(nansum(wgtTempReg_rad)))./squeeze(nansum(nansum(areaReg)));
    tempRegWgtAvg_full   = squeeze(nansum(nansum(wgtTempReg_full)))./squeeze(nansum(nansum(areaReg)));

    % -- Compute error bars based on standard error -- %
    
    if caseChoice==3
        clearvars numberPoints
        diffRainReg_std = nanstd(diffRain_phys(:)); 
        numberPoints    = find(isnan(diffRain_phys)==0);
        diffRainReg_stdErr = diffRainReg_std./sqrt(numel(numberPoints).*ratio_eff_full(ireg,1));
        
        clearvars numberPoints 
        diffTempReg_std = nanstd(diffTemp_phys(:));
        numberPoints    = find(isnan(diffTemp_phys)==0);
        diffTempReg_stdErr = diffTempReg_std./sqrt(numel(numberPoints).*ratio_eff_full(ireg,2)); 
        
        regionalMeans(ireg,:) =  [rainRegWgtAvg_phys, tempRegWgtAvg_phys];
        regionalErrorBars(ireg,:) = [diffRainReg_stdErr, diffTempReg_stdErr];
        
    elseif caseChoice==2
        clearvars numberPoints
        diffRainReg_std = nanstd(diffRain_rad(:)); 
        numberPoints    = find(isnan(diffRain_rad)==0);
        diffRainReg_stdErr = diffRainReg_std./sqrt(numel(numberPoints).*ratio_eff_full(ireg,1));
        
        clearvars numberPoints 
        diffTempReg_std = nanstd(diffTemp_rad(:));
        numberPoints    = find(isnan(diffTemp_rad)==0);
        diffTempReg_stdErr = diffTempReg_std./sqrt(numel(numberPoints).*ratio_eff_full(ireg,2)); 
        
        regionalMeans(ireg,:) =  [rainRegWgtAvg_rad, tempRegWgtAvg_rad];
        regionalErrorBars(ireg,:) = [diffRainReg_stdErr, diffTempReg_stdErr];
        
    elseif caseChoice==4 
        clearvars numberPoints
        diffRainReg_std = nanstd(diffRain_full(:)); 
        numberPoints    = find(isnan(diffRain_full)==0);
        diffRainReg_stdErr = diffRainReg_std./sqrt(numel(numberPoints).*ratio_eff_full(ireg,1));
        
        clearvars numberPoints 
        diffTempReg_std = nanstd(diffTemp_full(:));
        numberPoints    = find(isnan(diffTemp_full)==0);
        diffTempReg_stdErr = diffTempReg_std./sqrt(numel(numberPoints).*ratio_eff_full(ireg,2));   
        
        regionalMeans(ireg,:) =  [rainRegWgtAvg_full, tempRegWgtAvg_full];
        regionalErrorBars(ireg,:) = [diffRainReg_stdErr, diffTempReg_stdErr];
    end 
    
   
    clearvars ilon ilat areaReg imonthT imonthP rainPhys rainRad rainFull rainCtrl rain99phys rain99rad rain99full rain99ctrl
    clearvars diffRain_phys diffRain_rad diffRain_full tempPhys tempRad tempFull tempCtrl maskBoth_expand diffTemp_phys 
    clearvars diffTemp_rad diffTemp_full rain99phys_avg rain99rad_avg rain99full_avg rain99ctrl_avg wgtRainCtrl wgtRainPhys 
    clearvars wgtRainRad wgtRainFull rainCtrlRegWgtAvg rainPhysRegWgtAvg rainRadRegWgtAvg rainFullRegWgtAvg tempPhys_avg 
    clearvars tempRad_avg tempFull_avg tempCtrl_avg wgtTempCtrl wgtTempPhys wgtTempRad wgtTempFull tempCtrlRegWgtAvg tempPhysRegWgtAvg
    clearvars tempRadRegWgtAvg tempFullRegWgtAvg rainRegPhys rainRegRad rainRegFull tempRegPhys tempRegRad tempRegFull wgtRainReg_phys 
    clearvars wgtRainReg_rad wgtRainReg_full wgtTempReg_phys wgtTempReg_rad wgtTempReg_full rainRegWgtAvg_phys rainRegWgtAvg_rad 
    clearvars rainRegWgtAvg_full tempRegWgtAvg_phys tempRegWgtAvg_rad tempRegWgtAvg_full diffRainReg_std diffRainReg_stdErr 
    clearvars diffTempReg_std diffTempReg_stdErr 
    
  
    fprintf('Done with case %d \n', ireg);
end


                    
%Compute percent change of area-weighted average values
pctChangePhys = ((regionalPhys-regionalCtrl)./regionalCtrl).*100;
pctChangeRad  = ((regionalRad-regionalCtrl)./regionalCtrl).*100;
pctChangeFull = ((regionalFull-regionalCtrl)./regionalCtrl).*100;

%Plot regional bar charts
categories = categorical({'Rain','Soil Moisture','GPP',...
    'Net Sfc SW','Snowmelt','Runoff','ET'});

% -- Plot all bars in a single panel 

%figure;
newMeans = regionalMeans';
newErrors = regionalErrorBars'; 

% %Which percent changes are statistically significant? 
sig = NaN(size(newMeans));
sig(abs(newMeans)>=(2*newErrors))=1;
% 
% ctrs = 1:7;
% ax = axes;
% hBar = bar(ctrs,newMeans,1); 
% for k1=1:size(newMeans,2)
%    ctr(k1,:) = bsxfun(@plus, hBar(1).XData, [hBar(k1).XOffset]');
%    ydt(k1,:) = hBar(k1).YData;
% end
% hold on;
% errorbar(ctr',ydt',2.*newErrors,'k','linestyle','none','LineWidth',1.5)
% hold off
% legend('SE US','SE Asia','Central Australia')
% %legend('West Amazon','Central Africa','SE Australia')
% %legend('India','Horn of Africa')
% xticklabels(ax,{'Rain [mm/day]','Soil Moisture [kg/m^2]','GPP [g C/m^2/s]',...
%     'Net Sfc SW [10 W/m^2]','Snowmelt [mm/day]','Runoff [mm/day]','ET [mm/day]'});
% grid on
% set(gca,'FontSize',18)
% title('Multiply Stressed regions: PHYS-CTRL');
% ylabel('Absolute difference')
%  
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01, 0.1, 0.99, 0.9]);


%Save data that was used to create the bar plots above 
%save('/Users/meganfowler/gp_fuse/Flooding-physiology/MatlabData/CESMvarsControlling/variablesControlling_RADreg_PHYSexp.mat','regionalMeans','regionalErrorBars')

