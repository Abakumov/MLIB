function [dt, corr] = mycorr_hilbert(btrace, mtrace, Gdt)

% this function is designed to estimate timeshifts between traces
% it requires script fithyp.m
% Abakumov Ivan
% 4th February 2017

    if size(btrace,1)~=1
        btrace = btrace'; 
        mtrace = mtrace'; 
    end

    ub = btrace;
    um = mtrace;

    ub = ub - mean(ub);
    um = um - mean(um);
    
    ub = abs(hilbert(ub)); 
    um = abs(hilbert(um)); 
  
    J = xcorr(um,ub)/sqrt(sum(ub.^2)*sum(um.^2));
    
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



       



