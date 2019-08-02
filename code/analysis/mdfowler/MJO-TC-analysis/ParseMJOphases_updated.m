% Compute MJO phase from ERA-I observations (1985-2013) and split into
% favorable/unfavorable phases for a given ocean basin. 
% 
% Meg D. Fowler, 2017-05-02 
%


%% Determine MJO phase 
load('anomalies_ymean_1983-2013_ERA-I.mat');        %Loads array "anomalies_ymean"
load('timeForAnomalies_1983-2013_ERA-I.mat');       %Loads time array

[RMM1,RMM2,RMM1components,RMM2components] = calculate_RMMs_from_anomalies_WHeofs (anomalies_ymean);

[mjoindex,mjophase] = RMMstophasespace(RMM1,RMM2);

%Time from ERA-I RMMs 
timeVec = datevec(double(time));
yr      = timeVec(:,1);
mon     = timeVec(:,2);
day     = timeVec(:,3);

%% Read in Nino3.4 data and set to match RMM time period
% [ninoAnom, ninoAbs,ninoYr,ninoMon] = getNino34(); 
% 
% NinoTime = datenum(double(ninoYr),double(ninoMon),15*ones(size(ninoYr)));   %Monthly averages so day is set to mid-month arbitrarily
% 
% % Create "daily" array of Nino3.4 - each "day" of the month has the same
% % value, but the length will be the same as the RMM time vector 
% 
% ninoTimeDaily = NaN(size(time));
% ninoAnomDaily = NaN(size(time));
% ninoAbsDaily  = NaN(size(time)); 
% for i=1:numel(time)
%    vec   = squeeze(timeVec(i,:));   %Isolate time vector for day i
%    match = find(ninoYr==vec(1) & ninoMon==vec(2));  %Find where yr/mon of nino correspond to day i
%    
%    ninoTimeDaily(i) = NinoTime(match);  %Fill arrays for day i with yr/mon of nino data
%    ninoAnomDaily(i) = ninoAnom(match);
%    ninoAbsDaily(i)  = ninoAbs(match); 
% end


%% Define MJO events
%Set threshold for determining MJO events 
threshold = 1;

%Define MJO events as those above chosen threshold
iMJOevents    = find(mjoindex>=threshold); 
MJOeventIndex = mjoindex(iMJOevents);   %Extract indices of MJO events
MJOeventPhase = mjophase(iMJOevents);   %Extract phases of MJO events
MJOeventTime  = time(iMJOevents);       %Extract dates of MJO events

% MJOeventNinoAnom = ninoAnomDaily(iMJOevents); %Extract Nino3.4 conditions of MJO events
% MJOeventNinoAbs  = ninoAbsDaily(iMJOevents); 

MJOeventTimeVec = datevec(double(MJOeventTime));
MJOeventMonth   = MJOeventTimeVec(:,2);

%% Compare against WH observations... (debug option)
WHdebug=1;
if (WHdebug)
    %Read in text file with WH RMMs
    fid = fopen('rmmObs.74toRealtime.txt');
    data = textscan(fid, '%d %d %d %f %f %d %f %s\n', 'Headerlines',2,'Delimiter' , '	','EndOfLine','\r\n');
    fclose(fid);

    yr_obs    = data{1};
    mon_obs   = data{2};
    day_obs   = data{3};
    RMM1_obs  = data{4};
    RMM2_obs  = data{5};
    phase_obs = data{6};
    amp_obs   = data{7};
    timeObs   = datenum([double(yr_obs), double(mon_obs), double(day_obs)]);

    % Match with RMMs calculated in Mike's code... [longer time frame in WH]
    timeVec = datevec(double(time));
    yr      = timeVec(:,1);
    mon     = timeVec(:,2);
    day     = timeVec(:,3);

    iMatch    = find(yr_obs==yr(1) & mon_obs==mon(1) & day_obs==day(1));
    iMatchEnd = find(yr_obs==yr(end) & mon_obs==mon(end) & day_obs==day(end));
    timeObs   = timeObs(iMatch:iMatchEnd);
    %WH matched with our time frame
    RMM1_match = RMM1_obs(iMatch:iMatchEnd);
    RMM2_match = RMM2_obs(iMatch:iMatchEnd); 
    amp_match  = amp_obs(iMatch:iMatchEnd);
    phs_match  = phase_obs(iMatch:iMatchEnd); 

    %Read in OLR file to get lat/lon 
    OLRfile    = '/Volumes/MyPassport/Data/TCs/obs/NOAA_interpolated_OLR_1983-2013.nc';
    lat = ncread(OLRfile,'lat');
    lon = ncread(OLRfile,'lon');
    %Load coast file for plotting 
    coast = load('coast.mat');
    % Load anomalies of OLR and U for plotting (includes lat, which anomalies_ymean does not)
    load('OLRanom_new.mat');
    load('U200anomaly_withLat.mat');
    load('U850anomaly_withLat.mat');

    %Debug option to create diagnostic figure 
    event_debug = 0;
    if (event_debug)    
        [x,y] = meshgrid (lon,timeObs(it));
        it = find (timeObs >= datenum([2007 8 1 0 0 0]) & timeObs <= datenum ([2008 5 28 0 0 0]));  %MJO event

        figure;
        subplot (2,2,1); 
        %Match WH obs with our calculation of MJO index
        plot(timeObs(it),amp_match(it),'k')
        hold on;
        plot(time(it),mjoindex(it),'r')
        datetick('x','mmm-dd-yyyy','keepticks');
        datetick ('x','mm/dd','keepticks'); 
        legend('WH Obs','ERA-I')
        grid on;

        subplot (2,2,2); 
        %Hovmoller of OLR, with contours of U200 on top of that
        A = squeeze(nanmean(anomalyOLR(it,:,:),3)); % time x lon
        pcolor (x,y,A);
        shading flat; colorbar ('vert');  caxis ([-60 60]); 
        A = squeeze(nanmean(anomalyU200(it,:,:),3)); % time x lon
        hold on; 
        contour (x,y,A,[-8 -4 ],'k');   %Black contours for negative anom
        contour (x,y,A,[4 8],'w');      %White contours for postive anom
        datetick ('y','mm/dd','keepticks'); 

        subplot (2,2,3); 
        %Time evolution of MJO phase - expect steps, since phase can't change by
        %   more than one at a time. 
        %plot (RMM1(it),RMM2(it)); 
        plot (timeObs(it),mjophase(it));
        datetick ('x','mm/dd','keepticks'); xlabel ('Phase'); grid on;

        subplot (2,2,4);
        %RMM1 and RMM2 time evolution
        plot (timeObs(it),RMM1(it),'r'); hold on; plot (timeObs(it),RMM2(it),'b'); legend ({'RMM1','RMM2'}); 
        datetick ('x','mm/dd','keepticks'); 

        %error ('debug pause'); 
    end
end
%% Plot propogation of OLR (debug option) - requires above debug to be 1
prop_debug = 1;
if (prop_debug)    
    figure;    
    c = -33:0.5:33;
    for iphase=1:8
    %   it_mjo = find( mjoindex_winter>= threshold & mjophase_winter==iphase ); 
       it_mjo = find ((mon >= 6 | mon <= 11) & mjoindex' >= threshold & mjophase' == iphase); % Based on W-H RMMs.
       OLR_anom_map = squeeze(nanmean(anomalyOLR(it_mjo,:,:))); %Average of times matching condition
       subplot(8,1,iphase);
       contourf(lon,lat,OLR_anom_map',c,'LineColor','none');
       title(sprintf('June-Nov OLR anomaly: Phase %d',iphase));
       %title(sprintf('anomalyOLR: Phase %d in ours',iphase));
%        if iphase==1
%           title(sprintf('anomalyOLR: Phase %d in ours [threshold = %d, winter]',iphase,threshold));
%        end
       hold on; 
       set(gca,'FontSize',12);
       plot(coast.long, coast.lat, 'k','LineWidth',2);
       axis('xy','equal',[0 280 -20 20]);
       hcb=colorbar;
       caxis([-20 20])
       set(hcb,'FontSize',12)
       caxis([min(c) max(c)])
       
    end
           error ('debug pause'); 

end

%% Separate months and MJO phases for basin's TC season

basin    = {'ATL',   'NE_Pac', 'NW_Pac',  'N_Ind',     'S_Ind',     'S_Pac'}; 
TCseason = {[7:11],  [6:10],   [6:11],    [4:6,10:12], [1:4,11,12], [1:4,12]}; %Cells with months of the TC season of each basin 
favHurr   = {[2],    [1,7,8],  [6,7],     [4,5],       [3,4],       [7]};      %Favorable phases for hurricanes per basin
unfavHurr = {[7],    [5],      [3,4],     [7,8],       [7,8],       [2,3,5]};  %Unfavorable phases for hurricanes per basin
    
% Create a cell array that's month x MJO phase, so that each month of the 
%   TC season has separate instances of each MJO phase. 
for iMon = 1:12  
    for iPhase=1:8  
        imonth_phase = find(MJOeventMonth==iMon & MJOeventPhase'==iPhase); %Indices of where month and phase meet specification
        
        MJOindex_monthPhase{iMon,iPhase} = MJOeventIndex(imonth_phase);
        MJOtime_monthPhase{iMon,iPhase}  = MJOeventTime(imonth_phase); 
        NinoAnom_monthPhase{iMon,iPhase} = MJOeventNinoAnom(imonth_phase);
        NinoAbs_monthPhase{iMon, iPhase} = MJOeventNinoAbs(imonth_phase); 
    end
end

%% Compare Nino3.4 anomalies distributions in active phases for favorable/unfavorable MJO phase 

% NWPac_ninoFav      = NinoAnom_monthPhase(TCseason{3},favHurr{3});
% NWPac_ninoUnfav    = NinoAnom_monthPhase(TCseason{3},unfavHurr{3});
% NWPac_mjoFav       = MJOindex_monthPhase(TCseason{3},favHurr{3});
% NWPac_mjoUnfav     = MJOindex_monthPhase(TCseason{3},unfavHurr{3});
% NWPac_mjotimeFav   = MJOtime_monthPhase(TCseason{3},favHurr{3});
% NWPac_mjotimeUnfav = MJOtime_monthPhase(TCseason{3},unfavHurr{3});
% 
% SInd_ninoFav    = NinoAnom_monthPhase(TCseason{5},favHurr{5});
% SInd_ninoUnfav  = NinoAnom_monthPhase(TCseason{5},unfavHurr{5}); 
% SInd_mjoFav     = MJOindex_monthPhase(TCseason{5},favHurr{5});
% SInd_mjoUnfav   = MJOindex_monthPhase(TCseason{5},unfavHurr{5});
% SInd_mjotimeFav = MJOindex_monthPhase(TCseason{5},favHurr{5});
% SInd_mjotimeUnfav = MJOindex_monthPhase(TCseason{5},unfavHurr{5});
% 
% %Unpack cell arrays into a single long array for distribution checking
% NWPac_ninoFavFull   = [];
% NWPac_ninoUnfavFull = [];
% SInd_ninoFavFull    = [];
% SInd_ninoUnfavFull  = [];
% for i=1:numel(NWPac_ninoFav(:,1))    %Same number in fav/unfav and NWPac/SInd
%     for j=1:numel(NWPac_ninoFav(1,:))
%         NWPac_ninoFavFull   = [NWPac_ninoFavFull;   NWPac_ninoFav{i,j}];
%         NWPac_ninoUnfavFull = [NWPac_ninoUnfavFull; NWPac_ninoUnfav{i,j}];
%         
%         SInd_ninoFavFull    = [SInd_ninoFavFull;   SInd_ninoFav{i,j}];
%         SInd_ninoUnfavFull  = [SInd_ninoUnfavFull; SInd_ninoUnfav{i,j}];
%     end
% end
% 
%% Number of observations with strong ENSO events in fav/unfavorable 

% favNeg = find(SInd_ninoFavFull <= -1);
% favPos = find(SInd_ninoFavFull >= 1);
% unfavNeg = find(SInd_ninoUnfavFull <= -1);
% unfavPos = find(SInd_ninoUnfavFull >= 1);
% 
% % favNeg = find(NWPac_ninoFavFull <= -1);
% % favPos = find(NWPac_ninoFavFull >= 1);
% % unfavNeg = find(NWPac_ninoUnfavFull <= -1);
% % unfavPos = find(NWPac_ninoUnfavFull >= 1);

%% Create distributions

% nbins = 15;
% 
% figure;
% hold on;
% h1 = histogram(SInd_ninoFavFull,nbins,'Normalization','probability');
% h2 = histogram(SInd_ninoUnfavFull,h1.BinEdges,'Normalization','probability');
% title('SInd: Nino3.4 anomalies in MJO phases');
% legend('Favorable','Unfavorable')
% 
% 
%Fitting normal distribution
% pd_fav   = fitdist(NWPac_ninoFavFull,'Normal');
% pd_unfav = fitdist(NWPac_ninoUnfavFull,'Normal');
% xval     = -3:0.06:3;
% y_fav    = pdf(pd_fav,xval);
% y_unfav  = pdf(pd_unfav,xval);
% 
% figure; 
% plot(xval,y_fav,'b');
% hold on;
% plot(xval,y_unfav,'r');
% legend('Favorable phase','Unfavorable phase');
% title('NWPac: Kernel distribution of Nino3.4 anomalies in MJO phases');

%% Randomly undersample parts of the unfavorable phase in NWPac
% Nino3.4 distribution is uneven here. Strong El Nino's have 100 obs in
% unfavorable phase, but only 57 in the favorable one [for NWPac]. 

%%%%%------ Northwest Pacific -------------------------------------
% % Randomly remove where ENSO >=1 in all months... 
%  
% ratioFav = numel(favPos)/numel(favNeg); %Keep ratio of positive to negative events same in both fav/unfav phases
% numUnfPos = numel(unfavPos) - round(numel(unfavNeg).*ratioFav);     %Number of elements to remove
% 
% r = randperm(numel(unfavPos),numUnfPos); %all unique values w/ randperm
% 
% removal_indices = unfavPos(r);
% rKeep = 1:numel(NWPac_ninoUnfavFull);
% 
% toKeep = setdiff(rKeep,removal_indices);
% removedNWPac = NWPac_ninoUnfavFull(toKeep);
% 
% figure;
% hold on;
% h1 = histogram(NWPac_ninoFavFull,nbins,'Normalization','probability');
% h2 = histogram(removedNWPac,nbins,'Normalization','probability');
% title('NWPac: Histogram distribution of Nino3.4 anomalies in MJO phases...Random removal');
% legend('Favorable','Unfavorable');
% 
% figure;
% pd_fav   = fitdist(NWPac_ninoFavFull,'Normal');
% pd_unfav = fitdist(NWPac_ninoUnfavFull,'Normal');
% pd_unfavRem = fitdist(removedNWPac,'Normal');
% xval     = -3:0.06:3;
% y_fav    = pdf(pd_fav,xval);
% y_unfav  = pdf(pd_unfav,xval);
% y_unfavRem = pdf(pd_unfavRem,xval);
% 
% subplot(1,2,1);
% plot(xval,y_fav,'b');
% hold on;
% plot(xval,y_unfav,'r');
% legend('Favorable phase','Unfavorable phase');
% title('NWPac: Normal dist of Nino3.4anom');
% set(gca,'FontSize',18)
% 
% subplot(1,2,2);
% plot(xval,y_fav,'b');
% hold on;
% plot(xval,y_unfavRem,'r');
% legend('Favorable phase','Unfavorable phase');
% title('NWPac: Normal dist of Nino3.4anom...Random removal');
% set(gca,'FontSize',18)
% 
% set(figure(1),'Position',[97 285 820 590])

%%%%----------- South Indian -------------------------------------------
% Randomly remove where ENSO >=1 in favorable phase for all months... 
 
% ratioFav = numel(unfavPos)/numel(unfavNeg); %Keep ratio of positive to negative events same in both fav/unfav phases
% numUnfPos = numel(favPos) - round(numel(favNeg).*ratioFav);     %Number of elements to remove
% 
% r = randperm(numel(favPos),numUnfPos); %all unique values w/ randperm
% 
% removal_indices = favPos(r);
% rKeep = 1:numel(SInd_ninoFavFull);
% 
% toKeep = setdiff(rKeep,removal_indices);
% removedSInd = SInd_ninoFavFull(toKeep);
% 
% figure;
% hold on;
% h1 = histogram(removedSInd,nbins,'Normalization','probability');
% h2 = histogram(SInd_ninoUnfavFull,h1.BinEdges,'Normalization','probability');
% title('SInd: Nino3.4 anomalies...Random removal');
% legend('Favorable','Unfavorable');
% 
% figure;
% pd_fav   = fitdist(SInd_ninoFavFull,'Normal');
% pd_unfav = fitdist(SInd_ninoUnfavFull,'Normal');
% pd_unfavRem = fitdist(removedSInd,'Normal');
% xval     = -3:0.06:3;
% y_fav    = pdf(pd_fav,xval);
% y_unfav  = pdf(pd_unfav,xval);
% y_unfavRem = pdf(pd_unfavRem,xval);
% 
% subplot(1,2,1);
% plot(xval,y_fav,'b');
% hold on;
% plot(xval,y_unfav,'r');
% legend('Favorable phase','Unfavorable phase');
% title('NWPac: Normal dist of Nino3.4anom');
% set(gca,'FontSize',18)
% 
% subplot(1,2,2);
% plot(xval,y_fav,'b');
% hold on;
% plot(xval,y_unfavRem,'r');
% legend('Favorable phase','Unfavorable phase');
% title('NWPac: Normal dist of Nino3.4anom...Random removal');
% set(gca,'FontSize',18)
% 
% set(figure(1),'Position',[97 285 820 590])
% 
% 

%% Okay great, but now how do we actually remove those from the obs records? 

%%%----------------- Northwest Pacific -----------------------------------

% %Only removal is from NWPac_ninoUnfavFull/NWPac_mjoUnfav. 
% %   Based on removal_indices --> 42 observations need to be removed. 
% 
% %Array of indicies where each month/phase begin/end
% %   June/phs3, June/phs4, July/3, July/4, Aug/3, Aug/4, etc...
% istart = [1   81 166 211 264 312 347 373 431 461 520 599];
% iend   = [80 165 210 263 311 346 372 430 460 519 598 658];
% 
% %Unpacking the nino3.4 indices went in a different order than is used to
% %index the cell arrays [see notes, but it's unpacked one row at a time for
% %the nino indexing, and goes by column here for cell indexing]. 
% cell_ref = [1 7 2 8 3 9 4 10 5 11 6 12];
% 
% %Loop over each month/phase separation and find which removal indices fall
% %   within that range. 
% 
% new_NWPac_mjoUnfav{6,2} = NaN;
% new_NWPac_mjotimeUnfav{6,2}=NaN;
% for i=1:numel(istart)
%     irem1   = find(removal_indices>=istart(i) & removal_indices<=iend(i));
%     irem    = removal_indices(irem1);
%     unpack  = NWPac_mjoUnfav{cell_ref(i)};
%     unpack_time = NWPac_mjotimeUnfav{cell_ref(i)};
%     
%     dumKeep = 1:1:numel(unpack);
%     toKeep  = setdiff(dumKeep,irem-(istart(i)-1));
% 
%     kept    = unpack(toKeep);
%     kept_time = unpack_time(toKeep);
%     
%     new_NWPac_mjoUnfav{cell_ref(i)} = kept;  
%     new_NWPac_mjotimeUnfav{cell_ref(i)} = kept_time;
% end
% 
% %Save basin's MJO indices and dates during fav/unfav phases in TC season
% %save('NWPac_FavUnfav_mjoIndexTime','NWPac_mjoFav','NWPac_mjotimeFav','new_NWPac_mjoUnfav','new_NWPac_mjotimeUnfav');

%%%-------------- South Indian -------------------------------------------

%Only removal is from Sind_ninoFavFull/SInd_mjoFav. 
%   Based on removal_indices --> 42 observations need to be removed. 


% %Array of indicies where each month/phase begin/end
% %   Jan/phs3, Jan/phs4, Feb/3, Feb/4, Mar/3, Mar/4, etc...
% istart = [1   119 182 277 364 444 514 568 629 708 768 855];
% iend   = [118 181 276 363 443 513 567 628 707 767 854 922];
% 
% %Unpacking the nino3.4 indices went in a different order than is used to
% %index the cell arrays [see notes, but it's unpacked one row at a time for
% %the nino indexing, and goes by column here for cell indexing]. 
% cell_ref = [1 7 2 8 3 9 4 10 5 11 6 12];
% 
% %Loop over each month/phase separation and find which removal indices fall
% %   within that range. 
% 
% new_SInd_mjoFav{6,2} = NaN;
% new_SInd_mjotimeFav{6,2}=NaN;
% for i=1:numel(istart)
%     irem1   = find(removal_indices>=istart(i) & removal_indices<=iend(i));
%     irem    = removal_indices(irem1);
%     unpack  = SInd_mjoFav{cell_ref(i)};
%     unpack_time = SInd_mjotimeFav{cell_ref(i)};
%     
%     dumKeep = 1:1:numel(unpack);
%     toKeep  = setdiff(dumKeep,irem-(istart(i)-1));
% 
%     kept    = unpack(toKeep);
%     kept_time = unpack_time(toKeep);
%     
%     new_SInd_mjoFav{cell_ref(i)} = kept;  
%     new_SInd_mjotimeFav{cell_ref(i)} = kept_time;
% end

%Save basin's MJO indices and dates during fav/unfav phases in TC season
%save('SInd_FavUnfav_mjoIndexTime','new_SInd_mjoFav','new_SInd_mjotimeFav','SInd_mjoUnfav','SInd_mjotimeUnfav');

