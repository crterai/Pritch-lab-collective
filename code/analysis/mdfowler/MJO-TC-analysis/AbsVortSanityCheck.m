%Sanity check of vorticity using daily output
%
% Meg Fowler, 2017-03-30

testFile = '/Volumes/MyPassport/Data/TCs/spinup.cam.h1.1870-01-21-00000.nc';

lat       = ncread(testFile, 'lat');
lon       = ncread(testFile, 'lon');
lev       = ncread(testFile, 'lev');
ctrl_hyam      = ncread(testFile, 'hyam');
ctrl_hybm      = ncread(testFile, 'hybm');
time      = nctimenoleap(testFile);
ctrl_U    = ncread(testFile, 'U');           %[m/s]
ctrl_V    = ncread(testFile, 'V');           %[m/s]
ctrl_PS   = ncread(testFile, 'PS');          %[Pa]

P0        = 1000;                       %Reference pressure (hPa)
pnew      = [950,850,700,600,500,250,200];

fprintf('Now interpolating to set pressure levels \n');
ctrl_U1      = permute(ctrl_U, [4 3 2 1]);
ctrl_V1      = permute(ctrl_V, [4 3 2 1]);
ctrl_PS1     = permute(ctrl_PS, [3 2 1]);

ctrl_U_press = ncl4matlab('vinth2p',  ctrl_U1, ctrl_hyam, ctrl_hybm, pnew,    ctrl_PS1, 1, P0, 1, 'True');
ctrl_V_press = ncl4matlab('vinth2p',  ctrl_V1, ctrl_hyam, ctrl_hybm, pnew,    ctrl_PS1, 1, P0, 1, 'True');

%Compute relative vorticity 
fprintf('Now calculating relative vorticity \n');

ctrl_div = ncl4matlab('uv2dv_cfd', ctrl_U_press, ctrl_V_press,lat,lon,2);
ctrl_relVort = ncl4matlab('uv2vr_cfd', ctrl_U_press, ctrl_V_press,lat,lon,2);
%ctrl_relVort = ncl4matlab('uv2vrF', ctrl_U_press, ctrl_V_press);

ctrlVort_850 = squeeze(ctrl_relVort(:,2,:,:));

omega = 7.292e-5;   %Rotation of the earth [rad/s]

%Convert latitude to be in radians instead of degrees
latRad = lat.*(pi/180);
planetVort = 2*omega*sin(latRad);  % [rad/s]

absVortCtrl = NaN(size(ctrlVort_850));  %Define blank array for absolute vorticity values
for t=1:numel(time)
   for j=1:numel(lon)
       zetaCtrl = squeeze(ctrlVort_850(t, :, j));    %Relative vorticity
       
       absVortCtrl(t,:,j) = zetaCtrl + planetVort'; %Absolute vorticity
   end
end

%% Plots 
coast = load('coast.mat');

figure;
contourf(lon,lat,squeeze(ctrlVort_850(8,:,:)),-1e-4:1e-5:1e-4,'linecolor','none')
hold on;
plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k');
axis('xy','equal',[0 360 -45 45])
hold off;
title('Relative Vorticity');
colormap('jet');
colorbar('location','eastoutside')

figure;
contourf(lon,lat,squeeze(absVortCtrl(8,:,:)),-2e-4:1e-5:2e-4,'linecolor','none')
hold on;
plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k');
axis('xy','equal',[0 360 -45 45])
hold off;
title('Absolute Vorticity');
colormap('jet');
colorbar('location','eastoutside')
