function inside = is_inside_rotated_octahedron_v2(x, y, z, center, side_length, rot_matrix)
    
    % Find rotation matrix from the inverse of the MTEX rotation matrix
    %use the transpose of the orientation matrix, per MTEX requirements
    %https://mtex-toolbox.github.io/MTEXvsBungeConvention.html
    rotated_point = rot_matrix* [(center(1)-x); (center(2)-y); (center(3)-z)];

    % Check if the rotated point is within the rotated octahedron
    if abs(rotated_point(1)) + abs(rotated_point(2)) + abs(rotated_point(3)) <= side_length
        inside = true;
    else
        inside = false;
    end
end