function [] = make3DKirchhoff_mig(ACQ, G, project_folder, inicdp, ncdp)

% Abakumov Ivan
% January 2017
% University of Hamburg

    velmod = MLD([project_folder '/velmod.mat']); 
    image = zeros(size(velmod));
    filter = sin(linspace(0,pi/2,6)).^2; 
    Gold = oldGrid(G); 
            
    for i=1:ncdp
        cdp = inicdp + (i-1);
        datafile = [project_folder '/stack_cdp_' num2str(cdp) '.mat'];  
        if exist(datafile, 'file') ==2
            result = MLD(datafile); 
            U = ones(2000,1); 
            U(1:625) = result.stack;
            U(30:35) = U(30:35).*filter'; 

            S = [ACQ.CDPX(ACQ.CDP == cdp), ACQ.CDPY(ACQ.CDP == cdp), 0];  
            
            tti = FSM3D_down(Gold, S, velmod); 
            tti = 2*tti; 
            amp = 1./(tti+G.dt).^2; 
            tind = floor(tti/G.dt) + 1; 
    
            dtti = tti/G.dt - floor(tti/G.dt);    
            image = image + amp.*((1-dtti).*U(tind) + dtti.*U(tind+1));
        end
    end
    save([project_folder '/image_cdp_' num2str(inicdp) '.mat'], 'image', '-v7.3'); 
end
