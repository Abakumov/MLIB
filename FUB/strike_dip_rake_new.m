function p = strike_dip_rake_to_moment(strike,dip,rake)
% This function computes the 3x3 unit potency tensor p for a given set of strike,
% dip  and rake angles. This definition contradicts the definition of Aki
% and Richards, who defined the x-direction as North. 
% The definition used here uses x-direction as Easting, y-direction as
% Northing and z-Direction as Depth

% © Nepomuk Boitz, FU Berlin (boitz@geophysik.fu-berlin.de) July 2018


p(1,1) = sind(dip)*cosd(rake)*sind(2*strike) - sind(2*dip)*sind(rake)*cosd(strike)^2; 
p(2,2) = -(sind(dip)*cosd(rake)*sind(2*strike) + sind(2*dip)*sind(rake)*sind(strike)^2); 
p(3,3) = sind(2*dip)*sind(rake); 
p(1,2) = (sind(dip)*cosd(rake)*cosd(2*strike) + 0.5*sind(2*dip)*sind(rake)*sind(2*strike)); 
p(2,1) = p(1,2); 
p(1,3) = (cosd(dip)*cosd(rake)*sind(strike) - cosd(2*dip)*sind(rake)*cosd(strike)); 
p(3,1) = p(1,3); 
p(2,3) = (cosd(dip)*cosd(rake)*cosd(strike) + cosd(2*dip)*sind(rake)*sind(strike)); 
p(3,2) = p(2,3); 