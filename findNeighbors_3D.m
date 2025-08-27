function [neighborInd] = findNeighbors_3D(xInd,yInd,zInd,grid,topRow)
%FINDNEIGHBORS Function that returns the indices of the neighbors from the
%xInd, yInd position in grid

%Using a 28-neighbor basis (Moore Neighborhood)
neighborInd = [    -1    -1    -1;
    -1    -1     0;
    -1    -1     1;
    -1     0    -1;
    -1     0     0;
    -1     0     1;
    -1     1    -1;
    -1     1     0;
    -1     1     1;
     0    -1    -1;
     0    -1     0;
     0    -1     1;
     0     0    -1;
     0     0     1;
     0     1    -1;
     0     1     0;
     0     1     1;
     1    -1    -1;
     1    -1     0;
     1    -1     1;
     1     0    -1;
     1     0     0;
     1     0     1;
     1     1    -1;
     1     1     0;
     1     1     1];

%6-neighbor basis (Von Neumann neighborhood)
% neighborInd = [-1,0,0;1,0,0;
%      0,-1,0;0,1,0;
%      0,0,-1;0,0,1];

neighborInd=neighborInd+ [xInd,yInd,zInd];

origLength = length(neighborInd);

for i = 0:origLength-1
    
    %Reverse order indexing
    ind=origLength-i;
    
    %If the index is out of bounds, remove it
    if any(neighborInd(ind,:)==0) || any(neighborInd(ind,:)>length(grid)) || neighborInd(ind,3)>topRow
        neighborInd(ind,:)=[];
    end
end

%Convert from position to index
neighborInd = sub2ind(size(grid),neighborInd(:,1), neighborInd(:,2), neighborInd(:,3));

end

