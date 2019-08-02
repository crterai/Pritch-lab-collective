fileName= '/Users/meganfowler/gp_fuse/global_15min_physiologyControl/outflw1161.nc';

outflw = ncread(fileName,'outflw');
lat=ncread(fileName,'lat');
lon = ncread(fileName,'lon');

maskVal = zeros(size(outflw(:,:,1)));

for i=1:numel(lon)
    for j=1:numel(lat)
        if sum(isnan(outflw(i,j,:)))~=size(outflw(i,j,:))
            maskVal(i,j) = 1;
            %fprintf('Found valid gridcell!\n')
        end
    end
end

save('validGridMask_CaMa.mat','maskVal')