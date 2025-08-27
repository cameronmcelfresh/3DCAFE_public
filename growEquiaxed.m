function [grid,GrainCenters] = growEquiaxed(gridSize,numGrains)

%Construct a set of random grain centers
GrainCenters= rand(numGrains,3)*gridSize;

[xgrid,ygrid,zgrid]=meshgrid(1:gridSize,1:gridSize,1:gridSize);

xgrid=xgrid(:);
ygrid=ygrid(:);
zgrid=zgrid(:);


grid=zeros([gridSize,gridSize,gridSize]);

for i = 1:numel(grid)
    %calculate the distances to all grains
    dists = sqrt((xgrid(i)-GrainCenters(:,1)).^2 + ...
        (ygrid(i)-GrainCenters(:,2)).^2 + zgrid(i)-GrainCenters(:,3).^2);
    
    [~,closestGrain] = min(dists); %find the minimum distnace
    grid(i)=closestGrain;
end

grid=reshape(grid,[gridSize,gridSize,gridSize]);

end

