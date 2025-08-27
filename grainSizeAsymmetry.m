function [grainDetails] = grainSizeAsymmetry(gridSlice,dx,xgrid,ygrid)
%grainSizeAsymmetry Function to calculate the equivalent circle diameter
%and grain asymmetry. 

uniqueGrainID = unique(gridSlice);

grainDetails = zeros(length(uniqueGrainID),3);

%% Find the grain sizes
for i = 1:length(uniqueGrainID)
        grainDetails(i,1)=uniqueGrainID(i); %grain ID
        grainVol = sum(sum(gridSlice==uniqueGrainID(i)))*dx^2;
        grainDiameter = 2*sqrt(grainVol/3.1415); %equivalent circle diameter in meters
        grainDetails(i,2)=grainDiameter;
end
%% Find the principle axes
pAxes = findprincipleAxes(gridSlice,xgrid,ygrid,1,1);
pMax = max(pAxes(:,2:3),[],2);
pMin = min(pAxes(:,2:3),[],2);

grainDetails(:,3) = pMax./pMin;
grainDetails(grainDetails(:,3)>100,3)=1;

%remove any grain details with fewer than a given number of voxels
minVoxels=8;
for i = 1:length(uniqueGrainID)
        reverseCount = length(uniqueGrainID)-i+1;
        numVoxels = sum(sum(gridSlice==uniqueGrainID(reverseCount)));

        if numVoxels <minVoxels
            grainDetails(reverseCount,:)=[];
        end
end

end