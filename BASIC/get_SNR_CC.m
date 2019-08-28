function SNR = get_SNR_CC(X)


%% Author: Ivan Abakumov
% Freie Universität Berlin
% E-mail: abakumov_ivan@mail.ru
% Publication date: 26th August 2019

% size(X) = [length of trace x number of traces]

    [nt, ntrac] = size(X); 

    if nt < ntrac
        disp('Error: size of X should be: length of trace x number of traces')
    end

    C = zeros(ntrac);

    for i = 1:ntrac
        for j = 1:ntrac        
            [~,maxC] = mycorr(X(i,:), X(j,:), 1);
            C(i,j) = maxC; 
        end
    end

    nij = ntrac*(ntrac-1)/2;

    A = zeros(nij,ntrac);
    b = zeros(nij,1);

    k = 1; 
    for i=1:ntrac
        for j=(i+1):ntrac
            A(k,i) = 1; 
            A(k,j) = 1; 
            b(k) = -2*log(C(i,j));
            k = k+1; 
        end
    end

    S = inv(A'*A)*A'*b;  %#ok<MINV>
    Q = exp(S);

    SNR = nt./(Q-1);

end

