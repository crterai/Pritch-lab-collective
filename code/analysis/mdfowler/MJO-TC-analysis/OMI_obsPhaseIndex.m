%Read in daily OMI index and adapt Mike's code to compute MJO index and
%phase. 

function [mjoindex,mjophase,t2] = OMI_obsPhaseIndex()


filepath = '/Volumes/MyPassport/Data/OMI_EOFS/omi.1x.txt';

omidata = dlmread(filepath); %Read in data

yr    = omidata(:,1);
mon   = omidata(:,2);
day   = omidata(:,3);
pc1   = omidata(:,5);
pc2   = omidata(:,6);
pcSum = omidata(:,7); 

%% Code snippet from Mike's RMMstophasespace.m 
nphase = 8;

r2d = 180./(4.*atan(1.0));
ang = atan2(-pc1,pc2)*180/pi + 180;  % WARNING <--- this wasn't obvious!
% (but as above is only way to get the rmm1,rmm2 plot to match up with phase as 
% in GotSchalk)
mjoindex = sqrt(pc2.^2 + (-pc1).^2);
ang (ang < 0) = ang (ang< 0) + 360;
angbnd(1,:) = linspace (0,315,nphase);
angbnd(2,:) = linspace (45,360,nphase);
mjophase = NaN(size(pc1));
for iphase = 1:nphase
    it = find (ang >= angbnd(1,iphase) & ang <= angbnd(2,iphase));
    mjophase(it) = iphase;
end
mjoindex = sqrt(pc2.^2 + (-pc1).^2); 

% Check against Kiladis et al. figure for Dynamo period 
t = datetime(yr,mon,day);
t2 = datenum(t);
% start=datenum([2011 10 01 0 0 0]);
% finish=datenum([2011 10 31 0 0 0]);
% 
% idynamo = find(t2 >= start & t2<= finish);
% 
% dynIndex = mjoindex(idynamo);
% dynPhase = mjophase(idynamo);

end
