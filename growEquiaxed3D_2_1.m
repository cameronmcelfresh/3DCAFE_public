function [grid, GrainCenters] = growEquiaxed3D_2_1(gridSize, numGrains)
%% Function to create a starting grid from a voronoi tesselation - Much faster version!!
    % Construct a set of random grain centers
    GrainCenters = rand(numGrains, 3,'single') * gridSize; %convert to single precision

    [xgrid, ygrid, zgrid] = meshgrid(1:gridSize, 1:gridSize, 1:gridSize);
    positions = [xgrid(:), ygrid(:), zgrid(:)];
    positions = single(positions); %convert to single precision for memory

    % Calculate distances to all grains using vectorized operations
    %dists = sqrt(sum((permute(positions, [1, 3, 2]) - permute(GrainCenters, [3, 1, 2])).^2, 3));
    %[~, closestGrain] = min(dists, [], 2); % Find the index of the minimum distance

    dists = pdist2(positions, GrainCenters, 'euclidean');
    [~, closestGrain] = min(dists, [], 2); % Find the index of the minimum distance

    % Reshape the grid
    grid = reshape(closestGrain, [gridSize, gridSize, gridSize]);
end