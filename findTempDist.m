function [tempDist] = findTempDist(rowIter,sideLength,gridSize,beamPath,tempData)
%findTempDist Function to return the temperature distribution given the row
%of interest and the beam path

%Find the orientation vector 
if rowIter>1
    currentRow = rowIter;
    previousRow = rowIter-1;
else
    currentRow = rowIter+1;
    previousRow = rowIter;
end

orientVector = (beamPath(currentRow,2:3)-beamPath(previousRow,2:3));
orientVector(2)=orientVector(2)*-1; %flip y axis?
orientVector = orientVector/norm(orientVector); %normalize it
orientVector(1)=round(orientVector(1));


%Find the correct degree of rotation
rotDegree = [];
if all(orientVector==[1 0])
    rotDegree=0; %0 degrees
elseif all(orientVector==[0 1])
    rotDegree=1; %90 degrees CC
elseif all(orientVector==[-1 0])
    rotDegree=2; %180 degrees CC
elseif all(orientVector==[0 -1])
    rotDegree=3; %270 degrees CC
end

if isempty(rotDegree)
    orientVector
    rowIter
    pause(100);
end
tempDist = rot90(tempData.Temps,rotDegree);

%Find the matrix length needed to sample

heatFluxGridLength = tempData.dx*tempData.gridSize;%true size of the heat flux grid

conversionScale = sideLength/heatFluxGridLength;

lengthNeeded = round(tempData.gridSize*conversionScale/2); %find the % of grid needed to sample


%Amount to shift the beam by
beamX=round(beamPath(rowIter,2)/heatFluxGridLength*tempData.gridSize);
beamY=round(beamPath(rowIter,3)/heatFluxGridLength*tempData.gridSize);

%Shift the image - also include a shift back to the origin because the 
%beam spot has been centered in the grid
tempDist = circshift(tempDist,[beamY beamX]); %shift beam to true position
tempDist = circshift(tempDist,[-lengthNeeded -lengthNeeded]); %shift original position back to origin

%beamX
%beamY

%Sample the temp dist matrix
tempDist_sampled = tempDist(tempData.gridSize/2-lengthNeeded:tempData.gridSize/2+lengthNeeded,...
    tempData.gridSize/2-lengthNeeded:tempData.gridSize/2+lengthNeeded);

%Scale the image
tempDist = imresize(tempDist_sampled,[gridSize gridSize],'nearest');

end

