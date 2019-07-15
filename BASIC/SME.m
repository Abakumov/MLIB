function sdata=SME(data, okno)

% Smoth file
% Abakumov Ivan
% last update: 30 September 2014

mcase = length(size(data));


if mcase==2                 % 2D case
    [M,N] = size(data);
    sdata = zeros(M,N);
    wdata = ones(M,N);
    for i=1:M
        for j=1:N
            s = sum(sum(data(max(1, i-okno):min(M, i+okno), max(1, j-okno):min(N, j+okno))));
            w = sum(sum(wdata(max(1, i-okno):min(M, i+okno), max(1, j-okno):min(N, j+okno))));
            sdata(i, j) = s/w;
        end
    end
elseif mcase==3             % 3D case
    [L,M,N] = size(data);
    sdata = zeros(L,M,N);
    wdata = ones(L,M,N);
    for i=1:L
        for j=1:M
            for k=1:N
                s=sum(sum(sum( data(max(1, i-okno):min(L, i+okno), max(1, j-okno):min(M, j+okno), max(1, k-okno):min(N, k+okno) ))));
                w=sum(sum(sum(wdata(max(1, i-okno):min(L, i+okno), max(1, j-okno):min(M, j+okno), max(1, k-okno):min(N, k+okno) ))));
                sdata(i,j,k)=s/w;
            end
        end
    end
else
    disp('This code is designed only for 2D and 3D arrays');
end