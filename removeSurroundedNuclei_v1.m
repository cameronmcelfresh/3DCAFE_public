function [nuclei] = removeSurroundedNuclei_v1(nuclei,grid)
%removeSurroundedNuclei Function to remove nuclei that are surrounded

origLength = size(nuclei,1);

for nIter = 0:origLength-1

    %Reverse order indexing
    ind=origLength-nIter;

    %Remove the nuclei if the spot is being superheated by the laser
    if grid(nuclei(ind,4))<0
        nuclei(ind,:)=[];
        continue
    end

    [xC,yC,zC] = ind2sub(size(grid),nuclei(ind,4)); %find the cell the nuclei belongs to

    neighbors = findNeighbors_3D(xC,yC,zC,grid,length(grid)+1); %find the neighbor indices

   if all(grid(neighbors)>0) %remove the nuclei if it is surrounded
        nuclei(ind,:)=[];
   end

end