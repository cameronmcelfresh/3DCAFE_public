function [grid] = evolveGrid_3D(grid,topRow)
%evolveGrid_3D Function to grow the grid using CA in 3D

newGrid = grid; %make a grid to hold the changes that will occur during solidification

%Find all the points that are liquid
liqPoints = find(grid==0);

for index = 1:numel(liqPoints)
    
    %Find the location of the index
    [xInd,yInd,zInd] = ind2sub(size(grid),liqPoints(index));

    %Find the neighborhood of points
    [neighborInd] = findNeighbors_3D(xInd,yInd,zInd,grid,topRow);

    neighbors = grid(neighborInd); %find the grain ID of the neighbors

    %If surrounded by liquid, then skip the rest of the procedure
    if all(neighbors<1)
        continue;
    end

    neighbors = neighbors(neighbors>0); %ignore the liquid neighbors

    [M,~] = mode(neighbors); %find the mode and frequency

    %Make sure there aren't multiple modes
    if length(M)==1
        newGrid(liqPoints(index))=M;
    else
        newGrid(liqPoints(index))=M(0);
    end
end

grid = newGrid;

end