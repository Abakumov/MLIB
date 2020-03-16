function [R, theta] = get_Rotation_matrix(stress)
% Get rotation matrix Ry(theta)

    if (size(stress) == [3 3])
        theta = atan(2*stress(1,3)/(stress(1,1)-stress(3,3)))/2;
    elseif (length(stress) == 6)
        theta = atan(2*stress(5)/(stress(1)-stress(3)))/2;
    else
        disp('Error: wrong input')
    end

    R = zeros(3,3);
    R(1,1) = cos(theta);
    R(2,2) = 1; 
    R(3,3) = cos(theta);
    R(1,3) = +sin(theta);
    R(3,1) = -sin(theta);

end

