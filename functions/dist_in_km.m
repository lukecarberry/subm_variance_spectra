function [dx, dy] = dist_in_km(pt1, pt2)
% Distance between two points on the Earth, in nautical miles
% ptn must be in [lat, lon]

dphi = pt2(1) - pt1(1); % Latitude 
dtheta = pt2(2) - pt1(2); % Longitude in second position
phi = 0.5*(pt2(1) + pt1(1)); 

R = 6371; % km
dy = R*sin(dphi/180*pi);       
dx = R*cos(phi/180*pi)*sin(dtheta/180*pi); %km


end

