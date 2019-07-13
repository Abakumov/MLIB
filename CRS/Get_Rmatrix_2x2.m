function R = Get_Rmatrix_2x2(alpha, beta)
    % Abakumov Ivan
    % University of Hamburg
    % e-mail: abakumov_ivan@mail.ru

% Get rotation matrix from ray-centered to Cartisian coordinates
% alpha is a dip angle [ - 90; 90]
% beta is an azimuthal angle [-90 90]
% See equation 1.21
 % Abakumov, I. (2017). Systematic analysis of double-square-root-based stacking operators

PHI = [cos(beta),  -sin(beta); ...
       sin(beta),   cos(beta)];
   
THETA = [ cos(alpha), 0; ...
          0,          1];
     
R = PHI*THETA; 