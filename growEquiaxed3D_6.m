function [grid] = growEquiaxed3D_6(gridSize, numGrains)

%Construct a set of random grain centers
GrainCenters= rand(numGrains,3)*gridSize;

[xgrid,ygrid,zgrid]=meshgrid(1:gridSize,1:gridSize,1:gridSize);

xgrid=xgrid(:);
ygrid=ygrid(:);
zgrid=zgrid(:);

grid=zeros([gridSize,gridSize,gridSize]);

for i = 1:numel(grid)
    %calculate the distances to all grains
    dists = pdist2([xgrid(i),ygrid(i),zgrid(i)], GrainCenters);
    [~,closestGrain] = min(dists); %find the minimum distnace

    grid(i)=closestGrain;
end

grid=reshape(grid,[gridSize,gridSize,gridSize]);

end


