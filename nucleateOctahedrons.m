function [grid,nuclei] = nucleateOctahedrons(grid,nuclei,euler,t,topRow)
%nucleateOctahedrons Function to find the position of new octahedrons

%Find all the points that are liquid
liqPoints = find(grid==0);

for index = 1:numel(liqPoints)
    
    %Find the location of the index
    [xInd,yInd,zInd] = ind2sub(size(grid),liqPoints(index));

    %Find the neighborhood of points
    [neighborInd] = findNeighbors_3D(xInd,yInd,zInd,grid,topRow);

    neighbors = grid(neighborInd); %find the grain ID of the neighbors

    %If surrounded by liquid, skip the rest of the procedure
    if all(neighbors<1)
        continue;
    end

    if sum(neighbors>1)<10 %avoid corner cases of "edge" growth - remove
        continue;
    end

    if size(nuclei,1)>1 %Skip if a nuclei already exists on one of the neighbors
        if any(neighborInd'==nuclei(:,4))
            continue
        end
    end

    neighbors = neighbors(neighbors>0); %ignore the liquid neighbors

    [M,~] = mode(neighbors); %find the mode and frequency

    if rand>0.95
        %Make sure there aren't multiple modes
        if length(M)==1
                eulerIndex = find(euler(:,1)==M);
                grid(liqPoints(index))=M;
                nuclei = [nuclei;
                    xInd, yInd, zInd,...
                    liqPoints(index),...
                    euler(eulerIndex,2),euler(eulerIndex,3),euler(eulerIndex,4),...
                    t,...
                    0,...
                    M];
        else
                eulerIndex = find(euler(:,1)==M(0));
                grid(liqPoints(index))=M(0);
                nuclei = [nuclei;
                    XInd, yInd, zInd,...
                    liqPoints(index),...
                    euler(eulerIndex,2),euler(eulerIndex,3),nuclei(euler,4),...
                    t,...
                    0,...
                    M(0)];
        end
    end
end



end