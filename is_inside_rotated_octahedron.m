function inside = is_inside_rotated_octahedron(x, y, z, center, side_length, euler_angles)
    
    % Inverse rotation using Euler angles
    %R_x = [1, 0, 0; 0, cos(-euler_angles(1)), -sin(-euler_angles(1)); 0, sin(-euler_angles(1)), cos(-euler_angles(1))];
    %R_y = [cos(-euler_angles(2)), 0, sin(-euler_angles(2)); 0, 1, 0; -sin(-euler_angles(2)), 0, cos(-euler_angles(2))];
    %R_z = [cos(-euler_angles(3)), -sin(-euler_angles(3)), 0; sin(-euler_angles(3)), cos(-euler_angles(3)), 0; 0, 0, 1];
    %rotated_point = R_z * (R_y * (R_x * [(center(1)-x); (center(2)-y); (center(3)-z)]));
    
    rotmZYZ = eul2rotm(euler_angles,'ZYZ'); %use the built-in function instead
    rotated_point = rotmZYZ* [(center(1)-x); (center(2)-y); (center(3)-z)];

    % Check if the rotated point is within the rotated octahedron
    if abs(rotated_point(1)) + abs(rotated_point(2)) + abs(rotated_point(3)) <= side_length
        inside = true;
    else
        inside = false;
    end
end