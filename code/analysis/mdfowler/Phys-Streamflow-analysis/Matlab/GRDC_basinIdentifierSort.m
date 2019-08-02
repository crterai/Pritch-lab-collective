% Quickly check which river stations meet our criteria for validation
% against. 
%
% Megan D Fowler, 2016-09-19

[num txt raw] = xlsread('/Users/meganfowler/Documents/Irvine/Flooding/20160512_GRDC_Stations.xlsx','grdc_metadata');

raw       = raw(2:end,:);   %Remove header line
riverName = txt(2:end,6);

% grdcID        = num(:,1);
% lat           = num(:,9);
% lon           = num(:,10);
catchmentArea = num(:,11);
% dataYears     = num(:,17);
% startYear     = num(:,15);
% endYear       = num(:,16);

minArea  = 100000;  %Minimum catchment area required 

validData = find(catchmentArea>=minArea);

validBasins = unique(riverName(validData));

validRecords = raw(validData,:);

count=1;
for i=1:numel(validData)

   startDate = validRecords{i,15};
   endDate   = validRecords{i,16};
   
   rangeYears = startDate:endDate;
   
   k = find(rangeYears>=1970 & rangeYears<=2000);
   if numel(k)>=20       
      validAllData(count,:)   = validRecords(i,:);
       
      count = count+1;
   end

end
% uniqueBasin = find(unique(char(validBasins2)));
% names = char(validBasins2(uniqueBasin));
% 
% T=table(validBasins2(uniqueBasin)',validStartYear(uniqueBasin)',validEndYear(uniqueBasin)', ...
%     validArea(uniqueBasin)','VariableNames',{'River','Start','End','Area'});

filename='ValidBasins-newCriteria.txt';
%writetable(T,filename);

Tall = table(validAllData(:,1),validAllData(:,6),validAllData(:,8), validAllData(:,9),validAllData(:,10),...
    validAllData(:,11),validAllData(:,13),validAllData(:,15),validAllData(:,16),...
    validAllData(:,17),validAllData(:,18),'VariableNames',{'StationID','RiverName','CountryCode'...
    'lat','lon','catchmentArea','DownstreamStationID','YearStart','YearEnd','NumYears','PctMissing'});
writetable(Tall,filename);

validLons=cell2mat(validAllData(:,10));
validLats=cell2mat(validAllData(:,9));
validBasins=validAllData(:,6);

figure;
hold on;
coast = load('coast.mat');
plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k');
plot(validLons,validLats,'r*')
hold off;
title(['GRDC stations that fit new criteria (100,000km^2 area): ', num2str(numel(validLons))]);
axis([-180 180 -90 90]);




