function surrounded_indicator = isSurrounded(grainID,grid)
%ISSURROUNDED Function to tell whether a grain is completely surrounded by
%other grains or not 

surrounded_indicator=0; %0==not surrounded, 1==surrounded

neighborGrains = []; %collection of neighbor grains
grainInd = find(grid==grainID); % find the pixels in the grain
[r,~] = size(grainInd); %number of pixels in grain

grainSub = zeros([r,2]); %create a 2-column matrix to hold the indices
[grainSub(:,1),grainSub(:,2)] = ind2sub(size(grid),grainInd); %save the indices

%cycle through each of the pixels int he grain
for g = 1:r
    %collect all the neighbors
    pixelSub = findNeighbors(grainSub(g,1),grainSub(g,2),grid);
    
    %add the neighbor grains
    for p = 1:length(pixelSub)
        neighborGrains = [neighborGrains;grid(pixelSub(p,1),pixelSub(p,2))];
    end
end

%Check if all the neighbors are different grains
if all(neighborGrains>1)
    surrounded_indicator=1;
end

end

