function [grid] = growEquiaxed3D_3(gridSize, numGrains,layerHeight)
%% Function to create a starting grid from a voronoi tesselation - Much faster version!!
%layer height is the number of on each individual layer
%The total number of grains will be divided equally between each layer

numLayers = round(gridSize/layerHeight);
grainsPerLayer = round(numGrains/numLayers);

grid = zeros([gridSize gridSize gridSize]);

for layer = 1:numLayers
    % Construct a set of random grain centers
    GrainCenters = rand(grainsPerLayer, 3) * gridSize;
    GrainCenters(:,3) = rand(grainsPerLayer,1)*gridSize/numLayers;

    [xgrid, ygrid, zgrid] = meshgrid(1:gridSize, 1:gridSize, 1:layerHeight);
    positions = [xgrid(:), ygrid(:), zgrid(:)];

    % Calculate distances to all grains using vectorized operations
    dists = sqrt(sum((permute(positions, [1, 3, 2]) - permute(GrainCenters, [3, 1, 2])).^2, 3));

    [~, closestGrain] = min(dists, [], 2); % Find the index of the minimum distance

    % Reshape the grid
    nextLayer = reshape(closestGrain, [gridSize, gridSize, layerHeight]);

    %Add the next layer
    firstIter = 1+(layer-1)*layerHeight;
    lastIter = (layer)*layerHeight;
    grid(:,:,firstIter:lastIter)=nextLayer+grainsPerLayer*(1-layer);
end

%Randomly permute the grain ID throughout
allgrainIDs = unique(grid);
grainIDrand = randperm(length(allgrainIDs));
origGrid = grid;
for i = 1:length(allgrainIDs)
    grid(origGrid==allgrainIDs(i)) = grainIDrand(i);
end

end