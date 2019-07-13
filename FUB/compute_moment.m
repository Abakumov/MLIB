function M = compute_moment(p,C4d)
%% Computes Moment-Tensor from potency tensor and 4d-Stiffness tensor
% © Nepomuk Boitz, FU Berlin (boitz@geophysik.fu-berlin.de) January 2017

M = zeros(3); 
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
            M(i,j) = M(i,j) + C4d(i,j,k,l)*p(k,l); 
            end
        end
    end
end