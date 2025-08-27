function [grid,solidifyTime] = evolveGrid_3D_growthRate(grid,topRow,solidifyTime, dx, currentTime, growthRate)
%evolveGrid_3D_growthRate Function to grow the grid using CA in 3D, using
%the persribed growthRate (m/s)

newGrid = grid; %make a grid to hold the changes that will occur during solidification
newsolidifyTime = solidifyTime;

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

    %Find the neighbor wtih the most recent solidification time, which will
    %determine the minimum growth rate to cover the gap
    minGrowthRate = dx/(currentTime-min(solidifyTime(neighborInd)));

    if growthRate>minGrowthRate %only grow the cell if the growth rate is fast enough to cover the gap

        %Make sure there aren't multiple modes
        if length(M)==1
            newGrid(liqPoints(index))=M;
            newsolidifyTime(liqPoints(index))=currentTime; %set the current time to the solification time
        else
            newGrid(liqPoints(index))=M(0);
            newsolidifyTime(liqPoints(index))=currentTime; %set the current time to the solification time
        end
    else
        continue
    end

end

grid = newGrid;
solidifyTime= newsolidifyTime;


end