function velmod = voronoi2velmod(X,G,XX,ZZ)
 % Author: Abakumov Ivan
% Freie Universität Berlin
% E-mail: abakumov_ivan@mail.ru
% Publication date: 15th of May, 2018   
    VCx = X(1:40); 
    VCz = X(41:80); 
    VCv = X(81:120); 
        
    kdtree = KDTreeSearcher([VCx(:),VCz(:)]);
    f = @(x_in,y_in) VCv( kdtree.knnsearch([x_in(:),y_in(:)]) );
    
    velmod = f(XX(:),ZZ(:));
    velmod = reshape(velmod,G.nz,G.nx);
    velmod = velmod';



end

