function [d] = fit_grainIntersection(rotMat,grainZRotAng,a,b,dx)
%fit_grainIntersection Is a function that calculates the confinement
%directions for a grain given its:
%   orient: orientation matrix
%   s: slip directions
%   a,b: ellipse fitting parameters

%% Define the roation matrces
% xRot = [1 ,0, 0;
%     0 ,cosd(xRotAng),-sind(xRotAng);
%     0, sind(xRotAng), cosd(xRotAng)];
% 
% yRot = [cosd(yRotAng), 0 , sind(yRotAng);
%     0, 1, 0;
%     -sind(yRotAng), 0, cosd(yRotAng)];
% 
% zRot = [cosd(zRotAng), -sind(zRotAng), 0;
%     sind(zRotAng), cosd(zRotAng), 0;
%     0, 0, 1];

% rotMat = zRot*yRot*xRot; %Full rotation matrix for each grain

%% Define the grain rotation

grainRot = [cosd(grainZRotAng), -sind(grainZRotAng), 0;
    sind(grainZRotAng), cosd(grainZRotAng), 0;
    0, 0, 1];

%% Define the slip directions
%slip planes
%n = [1 1 1;
%    -1 1 1;
%    1 -1 1;
%    1 1 -1];

%burgers directions
s = [1 -1 0;
    0 1 1;
    1 0 1];

%% Calculate the intersected distance
%Array to hold the confinement distances
d = zeros([length(b),1]); 

bGrain = zeros(size(s));

%Change each slip direction to the coordinate axes of the grain
for i = 1:length(s)
    bGrain(i,:) =  s(i,:)*(rotMat);

    d(i) = 2/sqrt(abs(bGrain(i,1))/a^2 + abs(bGrain(i,2))/b^2);
end

end