function [Xk,Xn_app,U,S,mu,sigma] = get_pca(Xn,k)

% make PCA analysis

%% Author: Ivan Abakumov
% Freie Universität Berlin
% E-mail: abakumov_ivan@mail.ru
% Publication date: 2nd September 2019

% Input:
% Xn - data points in Rn space: Xn(i, :) = {x1, x2, x3, ..., xn}
% K - we will project the data from Rn to Rk, k < n 

% Output:
% Xk - data points in Rk space
% Xn_app - appxominate positions in Rn space
% U - matrix with eigenvectors
% S - matrix with eigenvalues

    [m, n] = size(Xn);
    
    if k>n
        disp('Error in PCA: k should be less than n')
    end

    % Step 1: Normalize the data
    mu = mean(Xn);
    X_norm = Xn - mu;

    sigma = std(X_norm);
    X_norm = X_norm./sigma;

    % Step 2: Compute covariation matrix and preform SVD

    Sigma = X_norm'*X_norm/m; 

    [U, S, ~] = svd(Sigma);

    % Step 3: project the data
    
    Xk = X_norm*U(:,1:k);

    % Step 4: find corresponding approximate points in Rn
    Xn_app = (Xk*U(:,1:k)').*sigma +mu; 
    
end
