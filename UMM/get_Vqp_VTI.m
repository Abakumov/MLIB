function Vqp = get_Vqp_VTI(Sample,delta)
% Calculate Vpq
% for exact formula see equation 4 and 6 in : 
% Tsvankin, I. (1996). P-wave signatures and notation for transversely 
% isotropic media: An overview. Geophysics, 61(2), 467-483.

    sinb2 = (sin(pi*Sample.Theta/180)).^2;
    cosb2 = (cos(pi*Sample.Theta/180)).^2;

    C11 = Sample.C11(1);
    C33 = Sample.C33(1);
    C44 = Sample.C44(1);
    rho = Sample.rho; 

    % see eq. 4
    C13 = sqrt(2*C33.*(C33-C44)*delta + (C33-C44).^2) - C44;

    % see eq. 6
    A = (C11 + C44).*sinb2 + (C33 + C44).*cosb2;
    B = (C11 - C44).*sinb2 - (C33 - C44).*cosb2;

    Vqp = sqrt(  (A + sqrt(B.^2 + 4*(C13 + C44).^2.*sinb2.*cosb2))/2/rho  );

end

