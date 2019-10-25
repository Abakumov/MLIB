function J  = costFunction_delta_angle(Sample,delta)
%COSTFUNCTION Compute cost function for least squares solution

    J = zeros(size(delta)); 
    
    for i = 1:length(delta)
        newSample = Sample;
        
        % part 1
        Theta = 0:1:360; 
        Psi = get_psi_VTI(Sample,delta(i),Theta);
        psiq = Sample.Theta; 
        thetaq = interp1(Psi,Theta,psiq,'pchip');
        newSample.Theta = thetaq; 
        % part 2
        ym = get_Vqp_VTI(newSample,delta(i));
        yo = Sample.Vqp;
        dy = ym - yo;    
        m = length(dy);
        J(i) = dy'*dy/m/2;
    end

end

