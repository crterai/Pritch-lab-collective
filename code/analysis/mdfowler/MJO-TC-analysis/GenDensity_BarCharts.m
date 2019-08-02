% Modify gendensity.m from Kerry's scripts to get average genesis density
% in five regions of interest for a given phase of the MJO. 
%
% Meg Fowler, 2019-05-02 

%% This is the part from Kerry's gendensity.m file.
%  (Hence the lack of comments)

params    %  Load parameters
%
clear latscat longscat y x z imin imax kmin kmax
clf('reset')
if exist('shearstore','var') == 0
    storm=[];
    shape='circ';
    load('temp.mat')
end   
load sorted
if strcmp(bas,'GB')
    projection=gproject;
end    
[nn,m]=size(vstore);
pifac=acos(-1)./180;
%
[~,jmax]=min(vstore,[],2);
jmax=jmax-1;
%
[~,jmin]=min(max((startv-vnet),0),[],2);
longscat=zeros(1,nn);
latscat=zeros(1,nn);
freqmask=zeros(1,nn);
if exist('yearstore','var') == 0
    yeartemp=2000+zeros(1,nn);
    freqtemp=freq;
else
    yeartemp=yearstore;
    freqtemp=freqyear;
    atemp=ismember((min(yearstore):max(yearstore)),yearstore);
    atemp=cast(atemp,'like',yearstore);
    yearset=nonzeros(atemp.*(min(yearstore):max(yearstore)));
end 
n=0;
for i=1:nn,
   if vmax(i) >= peakv
      n=n+1; 
      longscat(n)=longstore(i,jmin(i));
      latscat(n)=latstore(i,jmin(i));
      if (strcmp(bas,'MT')|| max(latstore(:,1)) > 50) && longscat(n) > 200
          longscat(n)=longscat(n)-360;
      end    
      [~,nyy]=ismember(yeartemp(i),yearset);
      freqmask(n)=freqtemp(nyy);
   end   
end
longscat(n+1:nn)=[];
latscat(n+1:nn)=[];
freqmask(n+1:nn)=[];
if strcmp(mapmode,'auto')
    xmin=min(nonzeros(longscat))-dellong;
    xmax=max(longscat)+dellong;
    ymin=min(nonzeros(latscat))-dellat;
    ymax=max(latscat)+dellat;
    xmin=gres*floor(xmin/gres);
    xmax=gres*ceil(xmax/gres);
    ymin=gres*floor(ymin/gres);
    ymax=gres*ceil(ymax/gres);    
else
    xmin=longmin;
    xmax=longmax;
    ymin=latmin;
    ymax=latmax;
end 
if strcmp(bas,'AL') && max(latstore(:,1)) < 55
    xmin=max(xmin,250);
    xmax=min(xmax,360);
    ymin=max(ymin,0);
%elseif strcmp(bas,'GB')
%    xmin=0;
%    xmax=360;
end 
x=xmin:gres:xmax;
y=ymin:gres:ymax;
an=max(size(y));
am=max(size(x));
z=zeros(an,am);
%
latr=gres*floor(latscat/gres);
latr=latr+mod(y(1),gres);
longr=gres*floor(longscat/gres);
longr=longr+mod(x(1),gres);
[latl,ay]=ismember(latr,y);
[longl,ax]=ismember(longr,x);
ax=max(ax,1);
ay=max(ay,1);
N=max(size(latl));
ga=freqmask.*latl.*longl./cos(pifac*latr);  % Corrected March, 2016
for i=1:N,
    z(ay(i),ax(i))=z(ay(i),ax(i))+ga(i);
end
z=z/(nn*gres^2);
zsyn=z;
%
% Calculate new map limits
%
imin=1;
imax=an;
kmin=1;
kmax=am;
if strcmp(mapmode,'auto')
    xmin1=xmax;
    xmax1=xmin;
    ymin1=ymax;
    ymax1=ymin;
    kmin=am;
    kmax=1;
    imin=an;
    imax=1;
    zcrit=0.01*max(max(z));
    for j=2:an-1;
        for k=2:am-1;
            if z(j,k) > zcrit
                xmin1=min(xmin1,x(k));
                xmax1=max(xmax1,x(k));
                ymin1=min(ymin1,y(j));
                ymax1=max(ymax1,y(j));
                kmin=min(kmin,k);
                kmax=max(kmax,k);
                imin=min(imin,j);
                imax=max(imax,j);
            end    
        end
    end  
    xmin=xmin1;
    ymin=ymin1;
    xmax=xmax1;
    ymax=ymax1;
end
%% TESTING: this is the part I'm trying to do - get regional averages 

%Get regional averages of each term 
% lonsLeft  = [100,120,130,140,160];
% lonsRight = [120,130,140,160,180];
% 
% ilats = find(y>=5 & y<=20);
% %Pick one of the five regions defined by lonsLeft and lonsRight above
% iReg=1;   
% ilons = find(x>=lonsLeft(iReg) & x<lonsRight(iReg)); 
% 
% z_Selection = z(ilats,ilons); %Isolate genesis density in that region (stored as z in Kerry's code)
% 
% %How can I eliminate all the land?? 
% landtemp = landcalc(x,y);     %Another of Kerry's scripts - defines land/ocean
% landMask = landtemp;
% landMask(landtemp==1) = NaN;  %Ignore land
% landMask(landtemp<1)  = 1;    %Keep ocean points 
% landMask_Selection = landMask(ilons,ilats); %Isolate region 
% 
% %Double check with a plot 
% figure;
% contourf(x(ilons),y(ilats),(z_Selection.*landMask_Selection'));
% hold on;
% plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k', 'LineWidth',2);
% axis([90 130 0 30])
% 
% maskedZ = z_Selection.*landMask_Selection'; %Get rid of the land part in genDensity
% 
% avg_Reg1 = nanmean(nanmean(maskedZ));

% COOL - so then all the steps that need to happen are: (1) isolate lon and
% lat region of interest; (2) mask out land; (3) take average over finite
% values 
%% Let's do it for real. 

%Regional definitions 
lonsLeft  = [100,120,130,140,160];
lonsRight = [120,130,140,160,180];
ilats = find(y>=5 & y<=20);

%Land Mask 
landtemp = landcalc(x,y);
landMask = landtemp;
landMask(landtemp==1) = NaN; 
landMask(landtemp<1)  = 1; 

%Loop over all five regions
for iReg=1:5 
    ilons = find(x>=lonsLeft(iReg) & x<lonsRight(iReg)); 
    z_Selection = z(ilats,ilons); 
    landMask_Selection = landMask(ilons,ilats); 
    maskedZ = z_Selection.*landMask_Selection';
    avg_Phs8(iReg) = nanmean(nanmean(maskedZ));
    
    clearvars ilons z_Selection landMask_Selection maskedZ 
end

%% Bar chart with just genesis density, similar to GPI one

% Instructions: 
%   Run above script with prep for each phase of the MJO in the MIT model downscaling. 
%   Clear all variables except the avg_PhsX (clearvars -except <var1> <var2>...). Then save all of those into a
%   larger array, "allPhases", as below. 

allPhases = [avg_Phs1; avg_Phs2; avg_Phs3; avg_Phs4; avg_Phs5; avg_Phs6; avg_Phs7; avg_Phs8]; %Save all individual phs genesis density into larger array
phaseMeans = mean(allPhases);
allPhases_norm = allPhases./phaseMeans;    %Normalize by mean genesis density in each region 

figure; 
regionLabels = {'100?-120?','120?-130?','130?-140?','140?-160?','160?-180?'};

%Plot bars and change colors
hBar = bar(1:8,allPhases_norm,'FaceColor','flat');
legend(regionLabels,'fontsize',12)
hBar(1).CData = [1 0 0];
hBar(2).CData = [0.9290 0.6940 0.1250];
hBar(3).CData = [0.4660 0.6740 0.1880];
hBar(4).CData = [0 0.4470 0.7410];
hBar(5).CData = [0.4940 0.1840 0.5560];
%Plot options
hold on;
title('MIT Model Genesis Density', 'fontsize',20); 
ylabel('Normalized Genesis Density', 'fontsize',16); 
xlabel('MJO Phase (OMI)', 'fontsize',16);
%ylim([0 0.055])
%Add lines to separate phases
plot([1.5 1.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')
plot([2.5 2.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')
plot([3.5 3.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')
plot([4.5 4.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')
plot([5.5 5.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')
plot([6.5 6.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')
plot([7.5 7.5], ylim,'Color',0.5*ones(1,3),'HandleVisibility','off')

