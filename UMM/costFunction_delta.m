function J  = costFunction_delta(Sample,delta)
%COSTFUNCTION Compute cost function for least squares solution

    J = zeros(size(delta)); 

    for i = 1:length(delta)
        ym = get_Vqp_VTI(Sample,delta(i));
        yo = Sample.Vqp;
        dy = ym - yo;    
        m = length(dy);
        J(i) = sum(dy.^2)/m/2;
    end

end
