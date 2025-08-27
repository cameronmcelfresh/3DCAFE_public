function [closestCorner] = find_closest_octahedronCorner(x, y, z, center, side_length, euler_angles)
%find_closest_octahedronCorner Function to find the closest corner fo the
%octahedron and return the normal direction

    %All principle axes
    princpAxes = [1,0,0;
                -1, 0 0;
                0,1,0;
                0,-1,0;
                0,0,1;
                0,0,-1];

    
    % Inverse rotation using Euler angles
    R_x = [1, 0, 0; 0, cos(-euler_angles(1)), -sin(-euler_angles(1)); 0, sin(-euler_angles(1)), cos(-euler_angles(1))];
    R_y = [cos(-euler_angles(2)), 0, sin(-euler_angles(2)); 0, 1, 0; -sin(-euler_angles(2)), 0, cos(-euler_angles(2))];
    R_z = [cos(-euler_angles(3)), -sin(-euler_angles(3)), 0; sin(-euler_angles(3)), cos(-euler_angles(3)), 0; 0, 0, 1];
   
    %rotated_corner = R_z * (R_y * (R_x * [x; y; z]));
    rotated_corners = R_z * (R_y * (R_x * princpAxes'*side_length));

    rotated_corners = rotated_corners';
    
    shiftedPoint=[x-center(1),y-center(2),z-center(3)]; %position of point of interest relative to the center of the octohedron

    dists = sqrt(sum((shiftedPoint-rotated_corners).^2,2)); %distances for the selected point to each corner

    %Return the position of the closest corner
    [~,ind] = min(dists);

    closestCorner = rotated_corners(ind,:)+center; %position of the closest corner


end