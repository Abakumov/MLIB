function T = Get_traveltime_2D_MF(M, H, model)

% Calculate MF traveltime approximation
% Abakumov Ivan
% 7st April 2016
% abakumov_ivan@mail.ru
% University of Hamburg

T = zeros(size(M));  

alpha = model(1);                % deg 
Rnip = model(2);                 % km 
Rn   = model(3);                 % km
v0   = model(4);                 % km/s 

t0 = 2*Rnip/v0; 

for i=1:size(M,1); 
    for j=1:size(M,2); 

        if (Rnip==0 || Rn==0)
            fmoveout=0;
        else
            delx=M(i,j)-H(i,j);
            dely=M(i,j)+H(i,j);
            s=delx+dely+2*sin(alpha)*delx*dely/Rnip;
            if (s==0.)
                s=0.01;
            end
            sigma=(delx-dely)/s;
            spl=Rnip+sigma*Rn;
            sml=Rnip-sigma*Rn;
            if (abs(spl)<0.0001)
                spl=sign(spl)*0.001;
            end
            if (abs(sml)<0.0001)
                sml=sign(sml)*0.001;
            end
            rp=Rn/spl*(1+sigma)*Rnip;
            rm=Rn/sml*(1-sigma)*Rnip;
            
            dlty=rm^2+rm*2*dely*sin(alpha)+dely^2;
            if (dlty<0.)
                dlty=0.;
            end
            delty=(sqrt(dlty)-rm)/v0;
            aaa=rm+dely*sin(alpha);
            if (aaa<0)
                delty=(-sqrt(dlty)-rm)/v0;
            end

            dltx=rp^2+rp*2*delx*sin(alpha)+delx^2;
            if (dltx<0.)
                dltx=0.;
            end
            deltx=(sqrt(dltx)-rp)/v0;
            bbb=rp+delx*sin(alpha);
            if (bbb<0.)
                deltx=(-sqrt(dltx)-rp)/v0;
            end
            fmoveout=deltx+delty;
        end
        T(i,j)=t0+fmoveout;
    end
end






























