function delta = get_delta_VTI(Sample,Vqp)
% Calculate Vpq
% for exact formula see equation 4 and 6 in : 
% Tsvankin, I. (1996). P-wave signatures and notation for transversely 
% isotropic media: An overview. Geophysics, 61(2), 467-483.

    C11 = Sample.C11(1);
    C33 = Sample.C33(1);
    C44 = Sample.C44(1);
    rho = Sample.rho; 

    theta = Sample.Theta/180*pi;
    sinb2 = (sin(theta)).^2;
    cosb2 = (cos(theta)).^2;

    A = (C11 + C44)*sinb2 + (C33 + C44)*cosb2;
    B = (C11 - C44)*sinb2 - (C33 - C44)*cosb2;
    % Vm
    C13pC442 = (  (2*rho*Vqp.^2 - A).^2 - B.^2)./(4*sinb2.*cosb2);  % (C11 + C44)^2

    delta = ( C13pC442 - (C33-C44).^2 )/(2*C33*(C33-C44));
    
end

