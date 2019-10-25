function Vqp = get_Vqp_VTI_weak_inclined(Sample,alpha,delta,epsilon,dTheta,mu)
% Calculate Vpq
% for exact formula see equation 16a in : 
% Thomsen, Weak anisotropy paper in Geophysics
 
    theta = (Sample.Theta + dTheta)/180*pi;
    theta = acos(cos(theta)*cos(mu/180*pi));
    sinb2 = (sin(theta)).^2;
    cosb2 = (cos(theta)).^2;

    Vqp = alpha*(1 + delta*sinb2.*cosb2 + epsilon*sinb2.^2);

end

