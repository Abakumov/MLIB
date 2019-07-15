G.nx=1000;
G.nz=1000;

%% PSO parameters
number_of_parameters=2;
number_of_particles=25;   
number_of_interations=160;   

c1=2.5586;
c2=1.3358; 
w =0.3925;

Xmin=diag([ 1,      1])*ones(number_of_parameters,number_of_particles); 
Xmax=diag([ G.nx,   G.nz])*ones(number_of_parameters,number_of_particles); 
Vmin=diag([-G.nx/3,-G.nz/3])*ones(number_of_parameters,number_of_particles); 
Vmax=diag([ G.nx/3, G.nz/3])*ones(number_of_parameters,number_of_particles); 

%% Initialization
X=Xmin+(Xmax-Xmin).*rand(number_of_parameters,number_of_particles); 
V=zeros(number_of_parameters,number_of_particles);

F=-((X(1,:)-350).^2 + (X(2,:)-550).^2); 

[gmax,maxpos]=max(F);
Xgmax=X(:,maxpos)*ones(1,number_of_particles);

pmax=-inf(1,number_of_particles);
Xpmax=X;

%% Iterations
for i=1:number_of_interations     
    V=w*V+c1*rand(number_of_parameters,number_of_particles).*(Xpmax-X)...
         +c2*rand(number_of_parameters,number_of_particles).*(Xgmax-X);
    V=min(max(V,Vmin),Vmax);
    X=X+V;  
    ind=((X<Xmin)+(X>Xmax))==1;
    Xrand=Xmin+(Xmax-Xmin).*rand(number_of_parameters,number_of_particles); 
    X(ind)=Xrand(ind);  
    F=-((X(1,:)-350).^2 + (X(2,:)-550).^2); 
    ind=F>pmax;
    pmax(ind)=F(ind);
    Xpmax(:,ind)=X(:,ind);
    [gmax,maxpos]=max(F);
    Xgmax=X(:,maxpos)*ones(1,number_of_particles);
end

%% Plot result
plot(X(1,:),X(2,:),'.k')