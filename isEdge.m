function is_edge_binary = isEdge(xIndex,yIndex,grid)
%isEdge Function to find potential heterogenous nucleation sites

is_edge_binary=0;

%ignore the simulation boundaries as heterogenous nucleation sites
if xIndex==1 || yIndex==1 || xIndex==length(grid) || yIndex==length(grid)
    return
end

if grid(yIndex,xIndex+1)==1
    is_edge_binary=1;
end

if grid(yIndex,xIndex-1)==1
    is_edge_binary=1;
end

% try
%     if grid(xIndex,yIndex+1)==1
%         is_edge_binary=1;
%     end
% end
% 
% try
%     if grid(xIndex,yIndex-1)==1
%         is_edge_binary=1;
%     end
% end

    
end

