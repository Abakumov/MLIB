function C4d = C2d_to_C4d(C2d)

C4d = zeros(3,3,3,3); 
for row=1:6 
    for col=1:6
        x = [1 2 3 3 3 2]; 
        y = [1 2 3 2 1 1]; 
        C4d(x(row),y(row),x(col),y(col)) = C2d(row,col); 
        C4d(x(row),y(row),y(col),x(col)) = C2d(row,col); 
        C4d(y(row),x(row),y(col),x(col)) = C2d(row,col); 
        C4d(y(row),x(row),x(col),y(col)) = C2d(row,col); 
    end
end
            
        