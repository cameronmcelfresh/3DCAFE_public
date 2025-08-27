function [growthIndices] = findNucleiGrowthSlice(xgrid,ygrid,growthRate,pos,dt,currentTime,nucleationTime)
%findNucleiGrowthSlice Function to find the indices of the "grid" that the
%nuclei will grow into during the dt iteration

startRadius = growthRate*(currentTime-nucleationTime-dt);
endRadius = growthRate*(currentTime-nucleationTime);

radialDist = sqrt((xgrid-pos(1)).^2 + (ygrid-pos(2)).^2);

growthMat = radialDist<endRadius & radialDist>startRadius;

growthIndices = find(growthMat==1);

end

