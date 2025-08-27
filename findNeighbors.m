function [neighborInd] = findNeighbors(xInd,yInd,grid)
%FINDNEIGHBORS Function that returns the indices of the neighbors from the
%xInd, yInd position in grid

%Using a 8-neighbor basis
neighborInd = [-1 0;-1 1;0 1;1 1;1 0;1 -1;0 -1;-1 -1] + [xInd,yInd];

origLength = length(neighborInd);

for i = 0:origLength-1
    
    %Reverse order indexing
    ind=origLength-i;
    
    %If the index is out of bounds, remove it
    if any(neighborInd(ind,:)==0) || any(neighborInd(ind,:)>length(grid))
        neighborInd(ind,:)=[];
    end
end

end

