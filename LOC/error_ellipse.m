function rellipse = error_ellipse(covariance)

% This function returns error ellipse corresponding to the covariance
% matrix. The function is based on the code of Vincent Spruyt, see the link
% 
% http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/
% 
% for the explanation how to draw the error ellipse for normally distributed
% data, given a chosen confidence value.
%
% The script usis distribution of chi-squared function. See matlab link for
% mode details. 
%
% https://de.mathworks.com/help/stats/chi2inv.html
%
% *Author*: Ivan Abakumov
%
% *Publication date*: 12th October 2017
%    
% *E-mail*: abakumov_ivan@mail.ru

%% Parameters

% Find a value that exceeds _confidence_interval_ of the samples from a 
% chi-square distribution with _number_of_deg degrees of freedom.
%
% Example: confidence_interval = 95%
%          number_of_defree = 2

confidence_interval = 0.688;
number_of_defree = 2;
chisquare_val = chi2inv(confidence_interval,number_of_defree);



%% Calculate the eigenvectors and eigenvalues
[eigenvec, eigenval] = eig(covariance);

% Get the largest eigenvector
leigenvec = eigenvec(:, 2);

% Get the largest and smallest eigenvalues
seigenval = eigenval(1,1);
leigenval = eigenval(2,2);

%% Calculate the angle between the x-axis and the largest eigenvector

phi = atan2(leigenvec(2), leigenvec(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(phi < 0)
    phi = phi + 2*pi;
end

theta_grid = linspace(0,2*pi);
a=sqrt(chisquare_val*leigenval);
b=sqrt(chisquare_val*seigenval);

% the ellipse in x and y coordinates 
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );

% Define a rotation matrix
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];

% Rotate the ellipse to angle phi
rellipse = [ellipse_x_r;ellipse_y_r]' * R;
