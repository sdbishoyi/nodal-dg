function [rhsPsi,rhsPi, rhsPhi] = ScalarWaveKerrRHS1D(Psi,Pi,Phi,rkerr,L,M,a,e,time)

Globals1D;

% Compute impedance
Zimp = ones(Np,K);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rH = 2*M;
%l = 0; the function has been tested for l = 0, 2

    V1 = (rkerr/rH).^(2*e);
    V2 = 1-(rH./rkerr).^(1+e);
    V3 = ((1+e*(1+ rH./rkerr))* L*(L+1))./(((1+e)^3).* rkerr.^2);
    V4 = (e+ (rH./rkerr)).^(1+e)./((rkerr.^2)*(1+e)^2); 
    pot = V1.*V2.*(V3+V4);     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define field differences at faces

dPi = zeros(Nfp*Nfaces,K); dPi(:) = Pi(vmapM)-Pi(vmapP);
dPhi = zeros(Nfp*Nfaces,K); dPhi(:) = Phi(vmapM)-Phi(vmapP);

Zimpm = zeros(Nfp*Nfaces,K); Zimpm(:) = Zimp(vmapM);
Zimpp = zeros(Nfp*Nfaces,K); Zimpp(:) = Zimp(vmapP);
Yimpm = zeros(Nfp*Nfaces,K); Yimpm(:) = 1./Zimpm(:);
Yimpp = zeros(Nfp*Nfaces,K); Yimpp(:) = 1./Zimpp(:); 


%%% Left boundary conditions
PiL= (Pi(1,1) -Phi(1,1))/2; dPi (1,1) = Pi(1,1) - PiL; 
PhiL =  (-Pi(1,1) +Phi(1,1))/2; dPhi (1,1) = Phi(1,1) - PhiL;


%%% Right boundary conditions
PiR= (Pi(end,end) + Phi(end,end))/2; dPi (end,end) = Pi(end,end) - PiR; 
PhiR =  (Pi(end,end) + Phi(end,end))/2; dPhi (end,end) = Phi(end,end) - PhiR;


% evaluate upwind fluxes
fluxPi = 1./(Zimpm + Zimpp).*(nx.*Zimpp.*dPhi - dPi);
fluxPhi = 1./(Yimpm + Yimpp).*(nx.*Yimpp.*dPi - dPhi);


% compute right hand sides of the PDEs
rhsPsi = -Pi;
rhsPi = (-rx.*(Dr*Phi) + LIFT*(Fscale.*fluxPi)) + pot.*Psi ;
rhsPhi = (-rx.*(Dr*Pi) + LIFT*(Fscale.*fluxPhi));

return

