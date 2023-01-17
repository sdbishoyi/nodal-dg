function [Psi_arr,Pi_arr,Phi_arr,t_arr] = ScalarWaveSch1D(Psi,Pi,Phi,L,M,a,e,FinalTime)

Globals1D;
time = 0;
%e = 0;
rH = 2*M;

%%%%%%%%%%%%%% Use Coordinate Conversion Routine here %%%%%%%%%%%%%%%%%%%%
if e ~=0
    rkerr = RstarToR_AOS(x,rH,e);
else
    rkerr = RstarToR(x,M);
%    rkerr = RstarToR_AOS(x,rH,0)
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Runge-Kutta residual storage  
resPi = zeros(Np,K); resPhi = zeros(Np,K); resPsi = zeros(Np,K); 

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));
CFL=0.2;  dt = CFL*xmin;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps;

Psi_arr={}; Psi_arr=[Psi_arr,Psi];
Pi_arr={}; Pi_arr=[Pi_arr,Pi];
Phi_arr={}; Phi_arr=[Phi_arr,Phi];
t_arr={}; t_arr=[t_arr,time];
% outer time step loop 

for tstep=1:Nsteps
   for INTRK = 1:5
      timelocal = time + rk4c(INTRK)*dt;
      [rhsPsi,rhsPi, rhsPhi] = ScalarWaveSchRHS1D(Psi,Pi,Phi,rkerr,L,M,a,e,timelocal);
       
      resPsi = rk4a(INTRK)*resPsi + dt*rhsPsi;     
      resPi = rk4a(INTRK)*resPi + dt*(rhsPi);
      resPhi = rk4a(INTRK)*resPhi + dt*(rhsPhi);

      Psi = Psi+rk4b(INTRK)*resPsi;     
      Pi = Pi+rk4b(INTRK)*resPi;
      Phi = Phi+rk4b(INTRK)*resPhi;
   
   end 
   % Increment time
   time = time+dt;
   Pi_arr=[Pi_arr,Pi];
   Phi_arr=[Phi_arr,Phi];
   Psi_arr=[Psi_arr,Psi];
   t_arr=[t_arr,time];
end
% time
% dt
return
