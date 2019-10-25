function J  = costFunction_delta_weak_inclined(Sample,alpha,delta,epsilon,dtheta,mu)
%COSTFUNCTION Compute cost function for least squares solution

    J = zeros(length(alpha),length(delta),length(epsilon),length(dtheta),length(mu)); 
    
    for i = 1:length(alpha)
        for j = 1:length(delta)
            for k = 1:length(epsilon)
                for l = 1:length(dtheta)
                    for n = 1:length(mu)
                        ym = get_Vqp_VTI_weak_inclined(Sample,alpha(i),delta(j),epsilon(k),dtheta(l),mu(n));
                        yo = Sample.Vqp;
                        dy = ym - yo;    
                        m = length(dy);
                        J(i,j,k,l,n) = sum(dy.^2)/m/2;
                    end
                end
            end
        end
    end
end
