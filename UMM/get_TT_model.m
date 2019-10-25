function TT = get_TT_model(G, Vp, Epsilon, Delta, Alfa, Tind, Source, Receiver)

TTI = 1e20*ones(G.nx, G.nz, length(Source.x)); 

for i = 1:25:length(Source.x)
    S = zeros(1,3);
    S(1) = Source.x(i);
    S(3) = Source.z(i);
    tti  = FSM2DTTI(G, S, Vp, Epsilon, Delta, Alfa);    
    tti(~Tind) = 1e20; 
    TTI(:,:,i) = tti;
end

TTI = TTI*1000;         % convert to micro sec
sTTI = TTI(:,:,151);    % single sensor
mTTI = min(TTI,[],3);   % multiple sensors 

TT.ssens = sTTI(Receiver.gx(1), Receiver.gz); 
TT.msens = mTTI(Receiver.gx(1), Receiver.gz); 






