function [dt, corr] = mywcorr(btrace, mtrace, fpeak, Gdt)

% this function is designed to estimate timeshifts between traces
% it requires script fithyp.m
% Abakumov Ivan
% 29th November 2016

    if size(btrace,1)~=1
        btrace = btrace'; 
        mtrace = mtrace'; 
    end

    ub = btrace;
    um = mtrace;

    ub = ub - mean(ub);
    um = um - mean(um);
   

    C = xcorr(um,ub)/sqrt(sum(ub.^2)*sum(um.^2));
    W = exp(-(linspace(-10,10,length(C))).^2); 
    L = (length(C)-1)/2; 
    tt = (-L:L)*Gdt; 
    w = 2*pi*fpeak; 
    C = C.*exp(-1i*w*tt); 
    J = abs(xcorr(C,W));

    [~, b] = max(J);
  
    if ( b==1 || b==length(J))
        dt = 0;
        corr = 0;
    else
        y1 = J(b-1);
        y2 = J(b);
        y3 = J(b+1);
        xmax = fithyp(b, y1, y2, y3);
       
        corr = J(b);
        dt = (xmax - (length(J)+1)/2)*Gdt;
    end
end



       



