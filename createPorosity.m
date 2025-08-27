function [gridPorous] = createPorosity(grid,xgrid,ygrid,zgrid,Tm,maxTemps,powderSize)
%createPorosity Function to create porosity based on the maximum
%temperature
%   Function returns grid with added LOF and entrapped gas porosity based
%   on the local maximum temperautres

gridPorous = grid;

% https://www.mathworks.com/matlabcentral/answers/405186-fill-area-with-random-circles-having-different-diameters
%% Create the LOF porosity
R_mean=powderSize;
R_STD=5;
maxIters=500000;
[C,R,~,~] = spherePacking(length(grid),R_mean,R_STD,maxIters); %create packed sphereical volume

subMeltInds = find(maxTemps<Tm);

for i=1:length(subMeltInds) %cycle through each point to see if it lies in a solid or void region
    ind = subMeltInds(i);
    dists = pdist2(C,[xgrid(ind),ygrid(ind),zgrid(ind)]); %distance from each point, in meters
    if all(dists>R) %convert radius from microns to meters
        gridPorous(ind)=0;
    end
end

%% Visualize Gas Porosity 
% probGasPorosity = exp(-(Tm-maxTemps)/Tm_alpha)*realGridSize/length(grid);
% porosityTrue = normrnd(mu_GP,STD_GP,size(grid)); %sample whether or not gas porosity occurs
% porositySize = normrnd(mu_GPSize,STD_GPSize,size(grid)); %size of porosity for each position;
% 
% %Cycle through the locations of new gas porosity and 
% newPorosityInds = find(probGasPorosity>porosityTrue);
% 
% for i = 1:length(newPorosityInds)
%     cellNum = newPorosityInds(i);
%     xLoc = xgrid(cellNum);
%     yLoc = ygrid(cellNum);
% 
%     gridPorous(sqrt((xgrid-xLoc).^2 + (ygrid-yLoc).^2)<porositySize(cellNum))=0;
% end
length(C)
title = {'x' 'y' 'z' 'radius'};
C = [title; num2cell(horzcat(C,R))];
writecell(C,"powder.csv");


end