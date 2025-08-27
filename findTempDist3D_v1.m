function [tempDist_final] = findTempDist3D_v1(rowIter,topRow,sideLength,gridSize,beamPath,tempData)
%findTempDist3D Function to return the temperature distribution given the row
%of interest and the beam path

beamX=beamPath(rowIter,2); %extract the x,y,z locations
beamY=beamPath(rowIter,3);
%beamZ=beamPath(rowIter,4);

beamUnitOrinetation=beamPath(rowIter,5:6); %unit vector for beam orientation
CosTheta = max(min(dot([1,0],beamUnitOrinetation)/(norm([1,0])*norm(beamUnitOrinetation)),1),-1);
rotDegree = real(acosd(CosTheta)); %rotation angle fo the beam

tempDist = imrotate_3D(tempData.Temps,rotDegree); % Rotate the temperature data, want to use imrotate3...

%Find the matrix length to sample
tempData.heatFluxGridLength = tempData.dx*tempData.gridSize;%true size of the heat flux grid
heatFluxGridLength = tempData.dx*tempData.gridSize;%true size of the heat flux grid
conversionScale = sideLength/heatFluxGridLength;
lengthNeeded = round(tempData.gridSize*conversionScale/2); %find the % of grid needed to sample - use +/- 1/2 of the length

%Amount to shift the beam by. Account for the FEA-simulated melt location
%as the center of the grid
beamX=round(beamX/heatFluxGridLength*tempData.gridSize-0.5*conversionScale*tempData.gridSize);
beamY=round(beamY/heatFluxGridLength*tempData.gridSize-0.5*conversionScale*tempData.gridSize);
%beamZ=round(beamZ/heatFluxGridLength*tempData.gridSize);

tempDist = circshift(tempDist,[beamY beamX 0]); %translate beam to true position

%Sample the temp dist matrix
tempDist_sampled = tempDist(round(tempData.gridSize/2)-lengthNeeded:round(tempData.gridSize/2)+lengthNeeded,...
    round(tempData.gridSize/2)-lengthNeeded:round(tempData.gridSize/2)+lengthNeeded,...
    end-2*lengthNeeded:end);

tempDist_sampled = imresizen(tempDist_sampled, gridSize/length(tempDist_sampled)); %resize the sampled temperatures to the appropriate length
tempDist_final = zeros([gridSize gridSize gridSize]); %grid to hold the temperatures
tempDist_final(:,:,1:topRow) = tempDist_sampled(:,:,end-topRow+1:end); %transfer the melt pool temperatures to the appropriate place in the grid

end