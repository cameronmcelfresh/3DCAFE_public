function [grid,nuclei] = findNewOctohedrons_v2(grid,nuclei,t,currentZHeight)
%findNewOctohedrons Function to the new octohedrons captured during growth

% column 1-3: x,y,z position of the nuclei
% column 4: index of the cell that the nuclei belongs to
% column 5-7: euler angles of the nuclei
% column 8: initial nucleation time of the nuclei
% column 9: size of the nuclei
% column 10: grainID of nuclei

% %Using a 28-neighbor basis (Moore Neighborhood)
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

numNuclei = size(nuclei,1);

for nIter = 1:numNuclei %loop through all nuclei and find the newly converted points

    [xSub,ySub,zSub] = ind2sub(size(grid),nuclei(nIter,4)); %find the xyz index of the nuclei
    possibleSites = neighborInd+[xSub,ySub,zSub]; %find the x,y,z col/row indicators 

    origLength = length(possibleSites);
    %Remove points that are outside of the grid
    for i = 0:length(possibleSites)-1
    
    %Reverse order indexing
    ind=origLength-i;
    
        %If the index is out of bounds, remove it
        if any(possibleSites(ind,:)==0) || any(possibleSites(ind,:)>length(grid)) || any(possibleSites(ind,3)>currentZHeight)
            possibleSites(ind,:)=[];
        end
    end

    %Check if any of the newly captured points are liquid
    possibleIndices = sub2ind(size(grid),possibleSites(:,1), possibleSites(:,2),possibleSites(:,3)); 
    %isLiquid = zeros(length(possibleIndices),1);

    isLiquid = grid(possibleIndices);

    %New points for octohedrons
    for i=1:length(isLiquid)

        if isLiquid(i)<1 %only create a new octahedron if the point is liquid

            %Check if the point is contained
            if ~is_inside_rotated_octahedron_v2(possibleSites(i,1), possibleSites(i,2), possibleSites(i,3), nuclei(nIter,1:3), nuclei(nIter,9), reshape(nuclei(nIter,11:end),[3,3])')
                continue
            end

            grid(possibleIndices(i))=nuclei(nIter,10); %convert the point on the grid

            %Find the closest corner for the new octahedron
            [closestCorner] = find_closest_octahedronCorner_v1(possibleSites(i,1), possibleSites(i,2), possibleSites(i,3), nuclei(nIter,1:3), nuclei(nIter,9), reshape(nuclei(nIter,11:end),[3,3])');

            %Center of the new nuclei, assume it has half the side length
            newCenter = (closestCorner-nuclei(nIter,1:3))/2+nuclei(nIter,1:3);

            %Add the new nuclei
            nuclei = [nuclei;...
                newCenter(1), newCenter(2), newCenter(3),... %x,y,z position
                possibleIndices(i),... %cell index that the nuclei belongs to
                nuclei(nIter,5),nuclei(nIter,6),nuclei(nIter,7),.... %euler angles of the nuclei
                t,... %nucleation time
                nuclei(nIter,9)/2,... %octahedron size
                nuclei(nIter,10),... %grain ID
                nuclei(nIter,11:end),... %rotation matrix in MTEX form
                ]; 
        end   
    end
end

end