% create_CaMa_netCDF
%
% This program takes raw daily runoff output from Gabe's simulations and
% sorts them into yearly files with the appropriate dimensions and
% specifications required for running through the CaMa-Flood model. This
% includes a bilinear cubic spline interpolation to get the model output to
% 0.5 degree, and switching longitudes from 0 to 360 to -180 to 180
% instead. 
%
% Megan D. Fowler, 2017-07-13
%   - Updated 2017-08-22 to read in a full 30 years, through an extra for
%     loop for year ranges. 
%
% !!! Must mount greenplanet first for the paths to be correct !!! 

runoffDir = '/Users/meganfowler/gp_fuse/Flooding-physiology/rawRunoff/';
cases = {'cesm1_0_6.1850_4xco2_fdbgb.1deg.002.clm2.h1.',... %4xCO2 atmosphere only 
         'cesm1_0_6.1850_4xco2_fixgb.1deg.001.clm2.h1.',... %4xCO2 land only
         'cesm1_0_6.1850_4xco2_fulgb.1deg.002.clm2.h1.',... %Full 
         'cesm1_0_6.1850_prei.1deg.001.clm2.h1.'};          %Control
% yrRange = {'0161-0170', '0171-0180',' 0181-0190'};
% yrRangeCtrl = {'0021-0030', '0031-0040',' 0041-0050'};
yrRange = {'0141-0150', '0151-0160'};
yrRangeCtrl = {'0001-0010', '0011-0020'};
tailStr = '.QRUNOFF.ncrcat.nc';

%CaseNames for saving output files
caseName = {'fdgb','fixgb','fulgb','prei'};
     
%Get land mask for zero-ing out ocean grid cells      
lmaskFile = [runoffDir,'cesm1_0_6.1850_4xco2_fulgb.1deg.002.clm2.h0.0181-0190.ncra.nc'];
landmask  = ncread(lmaskFile,'landmask');   %Land mask; 1=land, 0=ocean

%Read in default file used by CaMa for matching up with
%camaFile = '/Volumes/MyPassport/CaMa/runoff1990.nc';
camaFile = '/Users/meganfowler/gp_fuse/Flooding-physiology/MatlabData/runoff1990.nc';
dlon     = ncread(camaFile, 'lon');
dlat     = ncread(camaFile, 'lat'); 
dtime    = ncread(camaFile,'time');   %Literally a vector of days/year 1:365

for iCase = 4:4
    %if iCase==4
        %yrStarts = [1021, 1031, 1041];
    %else
        %yrStarts = [1161, 1171, 1181];
    %end
    
    %For first 20 years of downscaling:
    yrStarts = [0001,0011];   %Control 
    %yrStarts = [0141,0151];   %4xco2 experiments 
    
    for iYr = 1:2
        if iCase==4
            filename = [runoffDir, cases{iCase}, yrRangeCtrl{iYr}, tailStr]; 
        else
            filename = [runoffDir, cases{iCase}, yrRange{iYr}, tailStr];
        end

       runoff   = ncread(filename, 'QRUNOFF');  %[mm/s] 
       lon      = ncread(filename, 'lon'); 
       lat      = ncread(filename, 'lat');
       time     = nctimenoleap(filename);

       %Convert units on runoff to [mm/day]
       runoff = runoff*86400;

       %Fill ocean cells with 0's instead of NaN's (to match default input files
        %   used in CaMa-Flood model)
        for i=1:numel(lon)
            for j=1:numel(lat)
                if landmask(i,j)==0 && isnan(runoff(i,j,1))==1
                    runoff(i,j,:)=0;
                end
            end
        end

       %Handle negative runoff values by setting them to NaNs 
       runoff(runoff<0) = NaN;

       %Convert longitudes to be -180 to 180
       numlon = numel(lon);        %number of longitudes
       delt_lon = lon(2)-lon(1);  %longitude resolution in experiment file
       newlon = linspace(-180,(180-delt_lon),numlon)';   %Longitude array from -180 to just under 180

       ix_right = find (lon>180); 
       ix_left = find (lon <= 180);
       ix = [ix_right; ix_left]; 
       for iix = 1:length(ix)
            runoff2(iix,:,:) = runoff(ix(iix),:,:);
       end
       runoff = runoff2; 

       runoffFullRecord(:,:,:) = runoff;

       %%%---Interpolate to same grid as used in default file---%%%
       %%%---Split into yearly sections and create netCDF files--%%%
       nyrs = numel(time)/365;            %Number of years in experimental file
       daysperyear = 1:365;                %New time dimension 
       runoff = permute(runoff,[3 2 1]); %Rearrange so that dimensions are [time lev lat lon] or [time lat lon]

       start = 1;      %Start index for yearly indexing
       yearstr = yrStarts(iYr); %Start at some arbitrary year (1000+model-year)


       for t=1:nyrs
           yr_run  = runoff(start:365*t,:,:); %Isolate a single year of runoff

           new_run  = ncl4matlab('linint2',newlon,lat,yr_run,'True',dlon,dlat,0);    %Bilinear interpolation to new 1/2deg grid

           % Replace missing values with NaNs
           new_run(new_run   > 1E35) = NaN;

           %Rearrange so that dimensions are [lon lat time] again
           new_run = permute(new_run,[3 2 1]);    

           totalRecord(:,:,start:365*t) = new_run; 

           %Print progress message to screen
           fprintf('\n Completed re-gridding of data for year %s \n', num2str(t)); 

            %%%---Write out to netCDF file---%%%
            %Create netCDF file
            nc_filename = [caseName{iCase},'/runoff',num2str(yearstr),'.nc']; %Set filename 

            nccreate(nc_filename,'runoff',...       %Create variables
                      'Dimensions',{'lon',numel(dlon),'lat',numel(dlat),'time',365});
            nccreate(nc_filename,'lat','Dimensions',{'lat',numel(dlat)});
            nccreate(nc_filename,'lon','Dimensions',{'lon',numel(dlon)});
            nccreate(nc_filename,'time','Dimensions',{'time',365});

            %Write data to netCDF file
            ncwrite(nc_filename,'time',daysperyear);
            ncwriteatt(nc_filename,'time','long_name','time');
            ncwriteatt(nc_filename,'time','units','day');

            ncwrite(nc_filename,'lon',dlon);
            ncwriteatt(nc_filename,'lon','long_name','longitude');
            ncwriteatt(nc_filename,'lon','units','degrees_east');

            ncwrite(nc_filename,'lat',dlat);
            ncwriteatt(nc_filename,'lat','long_name','latitude');
            ncwriteatt(nc_filename,'lat','units','degrees_north');

            ncwrite(nc_filename,'runoff',new_run);  
            ncwriteatt(nc_filename,'runoff','long_name','total runoff (QRUNOFF)');
            ncwriteatt(nc_filename,'runoff','units','mm/day');
            ncwriteatt(nc_filename,'runoff','FillValue','1.e+20f');   

            %Print progress message to screen
            fprintf(' Created netCDF file for year %s \n \n', num2str(t)); 

            start = start+365;      %Increment start index by a year
            yearstr = yearstr+1;    %Increment year string by 1

       end
    end  
    fprintf(' Created netCDF files for case %s \n \n', cases{iCase});
end

% runoffYr = runoff(1:365,:,:);
% runoffYr_avg = squeeze(nanmean(runoffYr,1));
% newRun_Avg = squeeze(nanmean(new_run, 3));
%
% cRun = 0.05:0.1:25;
% coast=load('coast.mat');
% 
% figure;
% contourf(newlon,lat,squeeze(nanmean(squeeze(yr2raw),3))',cRun,'Linecolor','none')
% hold on;
% plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k', 'LineWidth',2);
% axis('xy','equal',[-180 180 -90 90]);
% title('Before interpolation');
% 
% figure;
% contourf(dlon,dlat,squeeze(nanmean(squeeze(yr2interp),3))',cRun,'linecolor','none')
% hold on;
% plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k', 'LineWidth',2);
% axis('xy','equal',[-180 180 -90 90]);
% title('After interpolation');

