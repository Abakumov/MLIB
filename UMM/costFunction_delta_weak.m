function J  = costFunction_delta_weak(Sample,alpha,delta,epsilon,dtheta)
%COSTFUNCTION Compute cost function for least squares solution

    J = zeros(length(alpha),length(delta),length(epsilon),length(dtheta)); 

    for i = 1:length(alpha)
        for j = 1:length(delta)
            for k = 1:length(epsilon)
                for l = 1:length(dtheta)
                    ym = get_Vqp_VTI_weak(Sample,alpha(i),delta(j),epsilon(k),dtheta(l));
                    yo = Sample.Vqp;
                    dy = ym - yo;    
                    m = length(dy);
                    J(i,j,k,l) = sum(dy.^2)/m/2;
                end
            end
        end
    end
end
