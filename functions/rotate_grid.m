function [latr,lonr,chlr] = rotate_grid(lat,lon,chl,rsublat,rsublon,dX)

[~, dyl] = dist_in_km([rsublat(1,1), rsublon(1,1)],[rsublat(1,2), rsublon(1,1)]);
[dxl, ~] = dist_in_km([mean(rsublat(1,:)), rsublon(1,1)],[mean(rsublat(1,:)), rsublon(1,2)]);
% [dxl, dyl] = dist_in_km([rsublat(1,1), rsublon(1,1)],[rsublat(end,2), rsublon(1,1)]);

clear latmrl lonmrl lonmr latmr
lat_lin = linspace(rsublat(1,1),rsublat(1,2),(dyl*1000)/dX);
lon_lin = linspace(rsublon(1,1),rsublon(1,2),(dxl*1000)/dX);

chl_obj = scatteredInterpolant(lon(:),lat(:),chl(:));
[lonr,latr] = meshgrid(lon_lin,lat_lin);

lonr = flipud(lonr);latr = flipud(latr);

chl_obj.Method = 'nearest';
chlr = chl_obj(lonr,latr);