% This script handles counting the total number of storms that fall into
% each hurricane category for a phase of the MJO, based on results from
% prep.m. 
% I am assuming that the variable "vnetmax" is indeed the one necessary to
% do this, based on comments in prep.m. I am also assuming the wind speeds
% are in knots, based on the units on potential intensity (vpstore) and
% units of returned wind speed in utransfull.m (pretty confident it's in
% knots). [Confirmed, it is]
%
% Meg Fowler, 2019-04-16 
% 

maxWind = max(vstore,[],2);
nanmean(maxWind)

cat1 = length(find(maxWind>=64 & maxWind<=82));
cat2 = length(find(maxWind>=83 & maxWind<=95));
cat3 = length(find(maxWind>=96 & maxWind<=112));
cat4 = length(find(maxWind>=113 & maxWind<=136));
cat5 = length(find(maxWind>=137));

fprintf(' Category 1: %i \n Category 2: %i \n Category 3: %i \n Category 4: %i \n Category 5: %i \n', cat1, cat2, cat3, cat4,cat5);

nbins = 63; 
edges = linspace(0, 250, nbins);
h=histogram(maxWind,edges);
title('PHASE 1: Max Lifetime Sfc Wind') 
hold on
plot([64 64], ylim,'k-','linewidth',2)
plot([83 83], ylim,'k-','linewidth',2)
plot([96 96], ylim,'k-','linewidth',2)
plot([113 113], ylim,'k-','linewidth',2)
plot([137 137], ylim,'k-','linewidth',2)

%   binCounts(1,:) = h.Values;
binCenters = h.BinEdges(1:(nbins-1)) / (h.BinWidth/2);

hold on; 
plot(binCenters,binCounts(1,:),'LineWidth',1.5);
plot(binCenters,binCounts(2,:),'LineWidth',1.5);
plot(binCenters,binCounts(3,:),'LineWidth',1.5);
plot(binCenters,binCounts(4,:),'LineWidth',1.5);
plot(binCenters,binCounts(5,:),'LineWidth',1.5);
plot(binCenters,binCounts(6,:),'LineWidth',1.5);
plot(binCenters,binCounts(7,:),'LineWidth',1.5);
plot(binCenters,binCounts(8,:),'LineWidth',1.5);
legend(['Phase 1'; 'Phase 2'; 'Phase 3'; 'Phase 4'; 'Phase 5'; 'Phase 6'; 'Phase 7'; 'Phase 8']);
xlabel('Average Storm Peak Windspeed');
ylabel('Histogram counts'); 

%% The above is super noisy; try instead fitting a Weibull dist to each array of max wind speeds

xTest = 0.5:.5:200; %400 points

maxWind = max(vstore,[],2);

parmhat = wblfit(maxWind);
pd = makedist('Weibull','a', parmhat(1), 'b', parmhat(2));
pdf_test = pdf(pd, xTest);

wbl_Phase(8,:) = pdf_test; 

figure;
hold on;
plot(xTest,wbl_Phase(1,:),'LineWidth',2)
plot(xTest,wbl_Phase(2,:),'LineWidth',2)
plot(xTest,wbl_Phase(3,:),'LineWidth',2)
plot(xTest,wbl_Phase(4,:),'LineWidth',2)
plot(xTest,wbl_Phase(5,:),'LineWidth',2)
plot(xTest,wbl_Phase(6,:),'LineWidth',2)
plot(xTest,wbl_Phase(7,:),'LineWidth',2)
plot(xTest,wbl_Phase(8,:),'LineWidth',2)
legend(['Phase 1'; 'Phase 2'; 'Phase 3'; 'Phase 4'; 'Phase 5'; 'Phase 6'; 'Phase 7'; 'Phase 8']);
xlabel('Storm peak wind speed','FontSize',14);
ylabel('Probability Density Function (PDF)','FontSize',14);
title('Weibull PDF of Max Storm Windspeed','FontSize',18);


