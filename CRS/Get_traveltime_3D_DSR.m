function t = Get_traveltime_3D_DSR(MM, HH, t0, w, M, N)
    % Abakumov Ivan
    % University of Hamburg
    % e-mail: abakumov_ivan@mail.ru
    % 3D DSR formula
    % see equation 1.19
    % size(w) = 2 x 1
    % size(M) = 2 x 2
    % size(N) = 2 x 2
    % size(MM) = 2 x fold
    % size(HH) = 2 x fold
    % size(t)  = 1 x fold
    %  ts = sqrt((t0+w'*(m-h))^2 + 2*t0*(m'*N*m - 2*m'*N*h + h'*M*h)); 
    %  tg = sqrt((t0+w'*(m+h))^2 + 2*t0*(m'*N*m + 2*m'*N*h + h'*M*h));

    
    fold = size(MM,2); 
    t = zeros(1, fold); 
    for i=1:fold
        m = MM(:,i); 
        h = HH(:,i); 
        mNm = m'*N*m; 
        mNh = 2*m'*N*h; 
        hMh = h'*M*h;
        ts = sqrt((t0+w'*(m-h))^2 + 2*t0*(mNm - mNh + hMh)); 
        tg = sqrt((t0+w'*(m+h))^2 + 2*t0*(mNm + mNh + hMh));
        t(i) = (ts + tg)/2; 
    end
end
    