function [saveF, saveX] = make_markov_chain(G,acq,XX,ZZ,texact,sigma,Xmin,Xmax,Vmin,Vmax,number_of_parameters,number_of_iterations)

% Abakumov Ivan

    saveF = zeros(number_of_iterations,1); 
    saveX = zeros(number_of_iterations,number_of_parameters); 

    X=Xmin+(Xmax-Xmin).*rand(number_of_parameters,1); 
    X(81:120) = 2000; 

    F = Misfit_bayesian_tomo(X,G,XX,ZZ,acq,texact,sigma);
 
    for i = 1:number_of_iterations     
        accept = 1; 
        while(accept)
            % make random update
            V=Vmin+(Vmax-Vmin).*rand(number_of_parameters,1); 

            % add the update (all if pdf is low, only one component if pdf is
            % high)
            if F > 0.1
                j = i - floor(i/120)*120 +1;
                Xn = X;
                Xn(j)=X(j)+V(j);  
            else
                Xn = X + V;
            end
            
            % find new pdf
            Fn = Misfit_bayesian_tomo(X,G,XX,ZZ,acq,texact,sigma);

            % accept this step or find a better one?
            P = Fn/F;
            if P>1
                accept=0;
                 F = Fn;
                 X = Xn; 
            else
                Prand = rand(1);
                if P>Prand
                    accept=0;
                    F = Fn; 
                    X = Xn; 
                end
            end
        end
        saveF(i,1) = F; 
        saveX(i,:) = X; 
    end


