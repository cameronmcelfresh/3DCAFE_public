function [tempDist_final] = findTempDist3D(rowIter,sideLength,gridSize,beamPath,tempData)
%findTempDist3D Function to return the temperature distribution given the row
%of interest and the beam path

beamX=beamPath(rowIter,2); %extract the x,y,z locations
beamY=beamPath(rowIter,3);
beamZ=beamPath(rowIter,4);

beamUnitOrinetation=beamPath(rowIter,5:6); %unit vector for beam orientation
CosTheta = max(min(dot([1,0],beamUnitOrinetation)/(norm([1,0])*norm(beamUnitOrinetation)),1),-1);
rotDegree = real(acosd(CosTheta)); %rotation angle fo the beam

tempDist = imrotate_3D(tempData.Temps,rotDegree); % Rotate the temperature data, want to use imrotate3...

%Find the matrix length to sample
tempData.heatFluxGridLength = tempData.dx*tempData.gridSize;%true size of the heat flux grid
heatFluxGridLength = tempData.dx*tempData.gridSize;%true size of the heat flux grid
conversionScale = sideLength/heatFluxGridLength;
lengthNeeded = round(tempData.gridSize*conversionScale/2); %find the % of grid needed to sample - use +/- 1/2 of the length

%Amount to shift the beam by
beamX=round(beamX/heatFluxGridLength*tempData.gridSize);
beamY=round(beamY/heatFluxGridLength*tempData.gridSize);
beamZ=round(beamZ/heatFluxGridLength*tempData.gridSize);

tempDist = circshift(tempDist,[beamY beamX beamZ]); %translate beam to true position

%Sample the temp dist matrix
tempDist_sampled = tempDist(round(tempData.gridSize/2)-lengthNeeded:round(tempData.gridSize/2)+lengthNeeded,...
    round(tempData.gridSize/2)-lengthNeeded:round(tempData.gridSize/2)+lengthNeeded,...
    round(tempData.gridSize/2)-lengthNeeded:round(tempData.gridSize/2)+lengthNeeded);

tempDist_final = imresizen(tempDist_sampled, gridSize/length(tempDist_sampled));

end