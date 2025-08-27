function [possibleIndices] = generateHetBoundaries(grid,beamPath,beamDiameter,realGridSize)
%GENERATEHETBOUNDARIES Function to generate the potential surfaces for
%heterogeneous growth.

possibleIndices=[];

%Find all the unique x values of the beam path
uniqueX = unique(beamPath(:,2));

% [X,~]=meshgrid((1:length(grid)),(1:length(grid)));
% X = X/length(grid)*realGridSize;

%All indices
allInd = reshape(1:length(grid)*length(grid),[length(grid),length(grid)]);

%Convert the beam diameter to pixel units
beamWidthPixels = beamDiameter/(2*realGridSize)*length(grid);

%Column iterator
xLoc = 1:length(grid);

%Loop through all the unique X values
for i = 1:length(uniqueX)
    
    %Find the location of the beam in pixel coordinates
    beamX=uniqueX(i)/realGridSize*length(grid);
    
    %Distance between the beam location and the x location
    dist = abs(beamX-xLoc);
    
    %Find all pixels directly outside of the beam
    beamColumns = find(dist>beamWidthPixels-8 & dist<(beamWidthPixels+12));
    %beamColumns = find(dist>beamWidthPixels-2 & dist<(beamWidthPixels+12));
    
    %Find the edge locations
    edgeLocs = allInd(:,beamColumns);
    
    %Add to the list
    possibleIndices=[possibleIndices,edgeLocs(:)'];
end

end

