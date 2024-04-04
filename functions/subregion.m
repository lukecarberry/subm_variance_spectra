function [latr,lonr,chlr] = subregion(lat,lon,chl,rsublat,rsublon)
% Luke Carberry, 03/16/2022
% lat, lon, chl all matrices that come from mat_ls_files.m

% sample sublat and sublon formats
% rsublat(1,:) = [33.55 35.2]; %Santa Barbara Channel 
% rsublon(1,:) = [-121.1 -118.9];

I = find(lat >= rsublat(1,1) & lat <= rsublat(1,2) & lon >= rsublon(1,1) & lon <= rsublon(1,2));
[y, x] = ind2sub(size(lon), I);
latr = lat(min(y):max(y),min(x):max(x));
lonr = lon(min(y):max(y),min(x):max(x));
chlr = chl(min(y):max(y),min(x):max(x));
