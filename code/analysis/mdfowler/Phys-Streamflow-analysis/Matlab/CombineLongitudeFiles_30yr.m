% PROJECT: Flooding Physiology
%
% Compile all longitude files into a single record for analysis. 
%
% Meg D. Fowler, 2017-09-19
%     Be sure to mount /gdata/pritchard2/mdfowler/ via sshfs


%% OPTIONS %%%
baseDir = '/Users/meganfowler/gp_fuse/Flooding-physiology/MatlabData/30yr/outflw_lonSeparated/'; 

caseName = 'control';   %Options: control, full, physiology, radiation
var      = 'outflw';    %Options: outflw, fldare 

%% Get list of files 
%fullDir = [baseDir,caseName,'/','fldare_avg','/'];   %Full directory for above options 
%fullDir = [baseDir,caseName,'/',var,'/'];   %Full directory for above options 
fullDir = baseDir %Set for case where need new outflw files with minYearly included 
cd( fullDir );
%files = dir('fldareAvg*.mat');
files = dir('*.mat');
fileList = {files.name};

%For fldareAvg, lon2 not saved - use lon2 from fldare
% fullDirLon = [baseDir,caseName,'/','fldare','/'];   %Full directory for above options 
% %fullDir = [baseDir,caseName,'/',var,'/'];   %Full directory for above options 
% cd( fullDirLon )
% filesLon = dir('*.mat');
% %files = dir('*.mat');
% fileListLon = {filesLon.name};

%% Load and save data from each file, with longitudes corresponding to those 
%  in a CaMa file. 

CaMaFile = '/Users/meganfowler/gp_fuse/Flooding-physiology/MatlabData/fldare1161.nc'; 

clon = ncread(CaMaFile,'lon');      %Longitudes in original data file
clat = ncread(CaMaFile, 'lat');     %Latitudes in original data file 
nyears = 30;                        %Number of years in analysis

%Define empty arrays to store data in 
perc90_full        = NaN(numel(clon),numel(clat));
perc99_full        = NaN(numel(clon),numel(clat));
perc999_full       = NaN(numel(clon),numel(clat));
avgVar_full        = NaN(numel(clon),numel(clat));
perc90yearly_full  = NaN(numel(clon),numel(clat),nyears);
perc99yearly_full  = NaN(numel(clon),numel(clat),nyears);
perc999yearly_full = NaN(numel(clon),numel(clat),nyears);
maxYearly_full     = NaN(numel(clon),numel(clat),nyears);
minYearly_full     = NaN(numel(clon),numel(clat),nyears); 
avgYearly_full     = NaN(numel(clon),numel(clat),nyears);

verbose = 1;
for iFile=1:numel(fileList)     %Loop over each file in directory
    fileName = fileList{iFile};
    %fileLon  = fileListLon{iFile};
    cd( fullDir )
    load(fileName);
    
%     cd( fullDirLon )
%     load(fileLon);
    
    fprintf('FldareAvg File: %s \n', fileName);
%     fprintf('FldareLon File: %s \n', fileLon);
    
    for j=1:numel(lon2)         %Loop over each lon in file
        ilon = find(clon==lon2(j)); %Match lon
        
        %Fill full arrays at matching longitude 
        perc90_full(ilon,:)  = perc90(j,:);
        perc99_full(ilon,:)  = perc99(j,:);
        perc999_full(ilon,:) = perc999(j,:);
        avgVar_full(ilon,:)  = avgVar(j,:);
        perc90yearly_full(ilon,:,:)  = perc90yearly(j,:,:);
        perc99yearly_full(ilon,:,:)  = perc99yearly(j,:,:);
        perc999yearly_full(ilon,:,:) = perc999yearly(j,:,:);
        maxYearly_full(ilon,:,:)     = maxYearly(j,:,:); 
        minYearly_full(ilon,:,:)     = minYearly(j,:,:); 
        avgYearly_full(ilon,:,:) = avgYearly(j,:,:);
    end
    
    if verbose==1
       fprintf('Assembled longitudes %f - %f \n', min(lon2),max(lon2));
    end
end

%% Save full files 
cd /Users/meganfowler/Documents/MATLAB/FloodProject/CaMaFlood/Physiology/30yr

%save([caseName,'_avgYearlyFldare','.mat'],'avgYearly_full');

save([caseName,'_',var,'_withMin.mat'],'perc90_full','perc99_full','perc999_full','avgVar_full',...
     'perc90yearly_full','perc99yearly_full','perc999yearly_full','maxYearly_full','minYearly_full','avgYearly_full');
