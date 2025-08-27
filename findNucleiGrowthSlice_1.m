function [growthIndices] = findNucleiGrowthSlice_1(xgrid,ygrid,growthRate,pos,dt,dx,grid,nucleiID)
%findNucleiGrowthSlice Function to find the indices of the "grid" that the
%nuclei will grow into during the dt iteration

nucleiX = xgrid(grid==nucleiID); 
nucleiY = ygrid(grid==nucleiID); %find the x-y grid positions that correspond to the nuclei

allRadialDists = sqrt((nucleiX-pos(1)).^2+(nucleiY-pos(2)).^2); 
startRadius = max(allRadialDists)*dx; %find the point furthest from the center, set it as the radius

growthStep = growthRate*dt; %distance grown during this dt
endRadius=startRadius+growthStep; %total radial distance grown

radialDist = sqrt((xgrid-pos(1)).^2 + (ygrid-pos(2)).^2)*dx;


growthMat = radialDist<endRadius;
growthIndices = find(growthMat==1);

end

