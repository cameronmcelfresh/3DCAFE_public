function inside = is_inside_rotated_octahedron_v1(x, y, z, center, side_length, euler_angles)
    
    % Find rotation matrix from the inverse of the MTEX rotation matrix
    %use the transpose of the orientation matrix, per MTEX requirements
    %https://mtex-toolbox.github.io/MTEXvsBungeConvention.html
    ori=orientation.byEuler(euler_angles);
    rotMat = ori.matrix';
    rotated_point = rotMat* [(center(1)-x); (center(2)-y); (center(3)-z)];

    % Check if the rotated point is within the rotated octahedron
    if abs(rotated_point(1)) + abs(rotated_point(2)) + abs(rotated_point(3)) <= side_length
        inside = true;
    else
        inside = false;
    end
end