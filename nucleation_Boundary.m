function isBoundary_binary = nucleation_Boundary(xIndex,yIndex,grid)
%isEdge Function to find potential heterogenous nucleation sites

if grid(yIndex,xIndex)~=0
    isBoundary_binary=0;
    return;
end

isBoundary_binary=0;

try
    if grid(yIndex,xIndex+1)>0
        isBoundary_binary=1;
        return;
    end
end

try
    if grid(yIndex,xIndex-1)>0
        isBoundary_binary=1;
        return;
    end
end

try
    if grid(yIndex+1,xIndex)>0
        isBoundary_binary=1;
        return;
    end
end

try
    if grid(yIndex-1,xIndex)>0
        isBoundary_binary=1;
        return;
    end
end

    
end

