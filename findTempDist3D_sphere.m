function [tempDist_final] = findTempDist3D_sphere(rowIter,sideLength,gridSize,beamPath,tempData)
%findTempDist3D_sphere Function to return the temperature distribution given the row
%of interest and the beam path. Assuming a spherical cutoff

beamRadius = 50e-6; %beam radius

beamX=beamPath(rowIter,2); %extract the x,y,z locations
beamY=beamPath(rowIter,3);
beamZ=beamPath(rowIter,4);

%Mesh grid positions and step size
[xgrid,ygrid,zgrid] = meshgrid((1:gridSize),(1:gridSize),(1:gridSize));
dx = sideLength/gridSize;

%Find the local distances
dist = sqrt( (xgrid*dx-beamX).^2 + (ygrid*dx-beamY).^2 + (zgrid*dx-beamZ).^2);

tempDist_final = ones([gridSize gridSize gridSize])*300; %set all the other temperatures to 300
tempDist_final(dist<beamRadius)=1500; %everything within the beam is 1500K

end

