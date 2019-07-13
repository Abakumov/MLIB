function [ t0, w, M, N ] = Get_model_stacking_parameters( X0, model )
%Get_model_stacking_parameters 
% This scrip estimates CRS parameters
% See equation 1.14 in thesis
% See web page http://en.wikipedia.org/wiki/Finite_difference
% for finite difference derivatives approximation
% Abakumov Ivan
% 27th August 2016
    
    dx = 0.1;
    dX = [dx 0 0; 0 dx 0; 0 0 dx]; 
        
    % t0
    t0 = Get_model_exact_traveltime(X0, X0, model); 
    
    % Vector w
    w = zeros(3, 1); 
    for i = 1:3
        w(i) = (Get_model_exact_traveltime( X0+dX(:,i),     X0+dX(:,i),     model )...  
              - Get_model_exact_traveltime( X0-dX(:,i),     X0-dX(:,i),     model ))/(2*dx);  
    end
        
    % Matrix M
    M = zeros(3); 
    for i = 1:3
        for j = 1:3
            M(i,j) = 0.5*(Get_model_exact_traveltime( X0-dX(:,i)-dX(:,j),     X0+dX(:,i)+dX(:,j),     model )...  
                        - Get_model_exact_traveltime( X0-dX(:,i)+dX(:,j),     X0+dX(:,i)-dX(:,j),     model )... 
                        - Get_model_exact_traveltime( X0+dX(:,i)-dX(:,j),     X0-dX(:,i)+dX(:,j),     model )... 
                        + Get_model_exact_traveltime( X0+dX(:,i)+dX(:,j),     X0-dX(:,i)-dX(:,j),     model ))/(4*dx*dx);  
        end
    end
        
    % Matrix N
    N = zeros(3); 
    for i = 1:3
        for j = 1:3
            N(i,j) = 0.5*(Get_model_exact_traveltime( X0+dX(:,i)+dX(:,j),     X0+dX(:,i)+dX(:,j),     model )...  
                        - Get_model_exact_traveltime( X0+dX(:,i)-dX(:,j),     X0+dX(:,i)-dX(:,j),     model )... 
                        - Get_model_exact_traveltime( X0-dX(:,i)+dX(:,j),     X0-dX(:,i)+dX(:,j),     model )... 
                        + Get_model_exact_traveltime( X0-dX(:,i)-dX(:,j),     X0-dX(:,i)-dX(:,j),     model ))/(4*dx*dx);  
        end
    end
end