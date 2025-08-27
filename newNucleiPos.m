function [centers,centerIndex] = newNucleiPos(grid,xgrid,ygrid,nucleationRate,nucleationRateSTD,dt,shapeFactor)
    %Function to return the positions of the new nuclei given the current state fo the simulation
    
    %Calc total number of new nuclei from homogenous nucleation
    newNuclei=round(normrnd(nucleationRate,nucleationRateSTD)*dt);
    
    %SETTING HOMOGENOUS NUCLEI TO ZERO
    %newNuclei=0;
    
    centers = [];
    centerIndex=[];
    
    openIndices = find(grid==0);
    
    if newNuclei>0 & length(openIndices)>newNuclei
        %Find a random permutation for homogenous nuclei
        primeIndex = randperm(length(openIndices),newNuclei);
        
        centerIndex = openIndices(primeIndex);
        
        centers = [xgrid(centerIndex),ygrid(centerIndex)];
        %centers = [xgrid(centerIndex);ygrid(centerIndex)]'; %old version? why?
        %centers = [primeIndex;primeIndex]';
    end
    
    %Update the new openIndices
    openIndices = setdiff(openIndices,centerIndex);
                    
    %Calc total number of new nuclei from heterogenous nucleation
    if shapeFactor<1
        
        centerIndexEdge=[];
        centersEdge=[];
        
        %Find the number of heterogeneous nuclei
        newHeteroNuclei=round(normrnd(nucleationRate/shapeFactor,nucleationRateSTD)*dt);
        %fprintf("%i hetero nuclei\n",newHeteroNuclei);
        
        %Specify all the edge locations
        edgeLocs = [];
        numEmptyEdges=0;
        
        %Find all the edge positions
        for kInd = 1:length(openIndices)
            xIndex = xgrid(openIndices(kInd));
            yIndex = ygrid(openIndices(kInd));
            
            if isEdge(xIndex,yIndex,grid)==1
                edgeLocs = [edgeLocs,openIndices(kInd)];
            end
        end
                
        if newHeteroNuclei>0 & length(edgeLocs)>newHeteroNuclei
            %Find a random permutation for homogenous nuclei
            primeIndex = randperm(length(edgeLocs),newHeteroNuclei);

            centerIndexEdge = edgeLocs(primeIndex);

            centersEdge = [xgrid(centerIndexEdge);ygrid(centerIndexEdge)]';
            %centers = [primeIndex;primeIndex]';
        end
           
        %Concatenate the new heterogenous nuclei
        if ~isempty(centerIndexEdge)
            centers = [centers;centersEdge];
            centerIndex = [centerIndex,centerIndexEdge];
        end
    end
            
end

