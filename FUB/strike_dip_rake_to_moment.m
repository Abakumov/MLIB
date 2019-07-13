function M = strike_dip_rake_to_moment(strike,dip,rake)
% This function computes the 3x3 Moment-tensor M for a given set of strike,
% dip  and rake angles. 

M(1,1) = - (sind(dip)*cosd(rake)*sind(2*strike) + sind(2*dip)*sind(rake)*sind(strike)^2); 
M(2,2) = sind(dip)*cosd(rake)*sind(2*strike) - sind(2*dip)*sind(rake)*cosd(strike)^2; 
M(3,3) = sind(2*dip)*sind(rake); 
M(1,2) = sind(dip)*cosd(rake)*cosd(2*strike) + 0.5*sind(2*dip)*sind(rake)*sind(2*strike); 
M(2,1) = M(1,2); 
M(1,3) = -(cosd(dip)*cosd(rake)*cosd(strike) + cosd(2*dip)*sind(rake)*sind(strike)); 
M(3,1) = M(1,3); 
M(2,3) = -(cosd(dip)*cosd(rake)*sind(strike) - cosd(2*dip)*sind(rake)*cosd(strike)); 
M(3,2) = M(2,3); 