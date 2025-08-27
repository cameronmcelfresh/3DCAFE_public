function [X,Y,Z,tempMatrix] = tempDistPreProcess(filename,centerPoint,unitOfLength,sizeToExtract,stepSize)
%tempDistPreProcess Extracts the temperature distribution from an arbitrary
%temperature distribtuion
%   filename: name of the file with the temperature distribution data (in
%    csv table format
%   centerPoint: A single vector [x,y,z] with the center point, to be
%     converted to 0,0,0
%   unitOfLength: the absolute length scale associated with the voxel
%     points, e.g. unitOfLength = 1e-6 then the x,y,z coordinates are in
%     microns
%   sizeToExtract: the outside bounds to stop sampling from the temperature
%      distribution

temps = readmatrix(filename);

temps(:,1)=[]; %remove the first column, which is just IDs

%Shift axes to new center
temps(:,1)=temps(:,1)-centerPoint(1);
temps(:,2)=temps(:,2)-centerPoint(2);
temps(:,3)=temps(:,3)-centerPoint(3);

[X,Y,Z] = meshgrid(-sizeToExtract/2:stepSize:sizeToExtract/2,...
    -sizeToExtract/2:stepSize:sizeToExtract/2,...
    -sizeToExtract:stepSize:0);

%Create an interpolation function
F = scatteredInterpolant(temps(:,1),temps(:,2)*1.5,temps(:,3)*1.25,temps(:,4));

% Interpolate the grid using the points
tempMatrix = F(X(:)./unitOfLength,Y(:)./unitOfLength,Z(:)./unitOfLength);

%Reshape into the appropriate matrix shape
tempMatrix = reshape(tempMatrix, size(X));

end