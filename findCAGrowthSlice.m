function [growthSlice,growthIndex] = findCAGrowthSlice(xIndex,yIndex,grid,gridSize,posInd)
%FINDCAGROWTHSLICE Function to find the growth slice using CA

growthSlice = [];
growthIndex = [];

%Check the +/-1 spots
possibleSteps = [-1 0;-1 1;0 1;1 1;1 0;1 -1;0 -1;-1 -1];

for step=1:length(possibleSteps)
    %If the grid is empty, fill it with the grain type of interest
    
    %Make sure we are not out of bounds
    if (yIndex+possibleSteps(step,2)>0 && yIndex+possibleSteps(step,2)<gridSize+1 && ...
            xIndex+possibleSteps(step,1)>0 && xIndex+possibleSteps(step,1)<gridSize+1)
            
        %IF the nearby pixels are not taken, take them
        if grid(yIndex+possibleSteps(step,2),xIndex+possibleSteps(step,1))==0
                    growthSlice = [growthSlice;
                        yIndex+possibleSteps(step,2),xIndex+possibleSteps(step,1)];
                    
                    growthIndex = [growthIndex;
                        posInd+possibleSteps(step,2)+possibleSteps(step,1)*gridSize];
        end
    end
end  


end

