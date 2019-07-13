function R = Get_Rmatrix_3x3(alpha, beta)

% Get rotation matrix from ray-centered to Cartisian coordinates
% alpha is a dip angle
% beta is an azimuthal angle

PHI = [cos(beta),  -sin(beta), 0; ...
       sin(beta),   cos(beta), 0; ...
       0            0          1];
   
THETA = [ cos(alpha), 0,  sin(alpha); ...
          0,          1,  0;          ...
         -sin(alpha), 0,  cos(alpha)];
     
R = PHI*THETA; 