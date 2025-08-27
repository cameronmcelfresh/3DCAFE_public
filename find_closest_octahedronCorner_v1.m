function [closestCorner] = find_closest_octahedronCorner_v1(x, y, z, center, side_length, rotMat)
%find_closest_octahedronCorner Function to find the closest corner fo the
%octahedron and return the normal direction

    %All principle axes
    princpAxes = [1,0,0;
                -1, 0 0;
                0,1,0;
                0,-1,0;
                0,0,1;
                0,0,-1];


    rotated_corners = rotMat * (princpAxes'*side_length);

    rotated_corners = rotated_corners';
    
    shiftedPoint=[x-center(1),y-center(2),z-center(3)]; %position of point of interest relative to the center of the octohedron

    dists = sqrt(sum((shiftedPoint-rotated_corners).^2,2)); %distances for the selected point to each corner

    %Return the position of the closest corner
    [~,ind] = min(dists);

    closestCorner = rotated_corners(ind,:)+center; %position of the closest corner


end