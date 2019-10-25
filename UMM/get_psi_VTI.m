function psi  = get_psi_VTI(Sample,delta,theta)

    C11 = Sample.C11(1);
    C33 = Sample.C33(1);
    C44 = Sample.C44(1);
    C13 = sqrt(2*C33.*(C33-C44)*delta + (C33-C44).^2) - C44;
    rho = Sample.rho; 

    theta = theta/180*pi;

    dtheta = 0.001; 

    % Vm
    sinb2 = (sin(theta-dtheta)).^2;
    cosb2 = (cos(theta-dtheta)).^2;
    A = (C11 + C44)*sinb2 + (C33 + C44)*cosb2;
    B = (C11 - C44)*sinb2 - (C33 - C44)*cosb2;
    Vm = sqrt(  (A + sqrt(B.^2 + 4*(C13 + C44)^2.*sinb2.*cosb2))/2/rho  );

    % Vp
    sinb2 = (sin(theta+dtheta)).^2;
    cosb2 = (cos(theta+dtheta)).^2;
    A = (C11 + C44)*sinb2 + (C33 + C44)*cosb2;
    B = (C11 - C44)*sinb2 - (C33 - C44)*cosb2;
    Vp = sqrt(  (A + sqrt(B.^2 + 4*(C13 + C44)^2.*sinb2.*cosb2))/2/rho  );

    V  = (Vp+Vm)/2; 
    dV = (Vp-Vm)/2/dtheta;

    
    

    psi = atan(dV./V)+theta;
    psi = psi*180/pi; 
end
