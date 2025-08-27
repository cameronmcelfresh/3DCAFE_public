function [nuclei] = removeSurroundedNuclei(nuclei,grid)
%removeSurroundedNuclei Function to remove nuclei that are surrounded

origLength = size(nuclei,1);

for nIter = 0:origLength-1

    %Reverse order indexing
    ind=origLength-nIter;

    [xC,yC,zC] = ind2sub(size(grid),nuclei(ind,4)); %find the cell the nuclei belongs to

    neighbors = findNeighbors_3D(xC,yC,zC,grid,length(grid)+1); %find the neighbor indices

   if all(grid(neighbors)>0) %remove the nuclei if it is surrounded
        nuclei(ind,:)=[];
    end

end