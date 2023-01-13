% function [psiarr,tarr]=ScalarWaveKerrDriver(L)
% Driver script for solving the 1D Maxwell's equations
Globals1D;

% Polynomial order used for approximation 
N=8; subd=240;
Np = N+1;
% Generate simple mesh

[Nv, VX, K, EToV] = MeshGen1D(-25,725,subd);
% Initialize solver and construct grid and metric
StartUp1D;

%%%%%%%%%%%%%%%%%%%% Solve Problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FinalTime = 750;
L=1;
M=0.5; a=0;
e = [0.1*(1 + 1e-4),0];
%e = 0.2*(1 + 1e-4);
%e=[0,0.1];

%%%%%%%%%%%%%% Initial conditions from Paper %%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigma=4; muu= 2; A = 20;
psi_in = A*exp(-(x-muu).^2/(2*sigma^2));
phi_in = A*exp(-(x-muu).^2/(2*sigma^2)).*(-(x-muu)./(sigma^2));
%pi_in = phi_in;
pi_in  = zeros(Np,K);    %(pi = 0 would be the initially static case.)

%%%%%%%%%%%%%%%% Gaurav's Initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sigma=4; muu= 2; A = 20;
% psi_in = zeros(Np,K);
% pi_in = A*exp(-(x-muu).^2/(2*sigma^2));
% phi_in = zeros(Np,K);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[psiarr_1,piarr_1,phiarr_1,tarr_1] = ScalarWaveKerr1D(psi_in,pi_in,phi_in,L,M,a,e(1),FinalTime);
tarr_1=cell2mat(tarr_1);

 psi_extr_1={};
 extr_posn = 250;
 [r,c]=find(abs(x-extr_posn)<=1);
 for i=1:length(tarr_1)
    psi_t=cell2mat(psiarr_1(i));
    psi=psi_t(r(1),c(1));
    psi_extr_1=[psi_extr_1,psi];
 end
 psi_extr_1= cell2mat(psi_extr_1);

 plot(tarr_1(4:end),log10(abs(psi_extr_1(4:end))),"LineWidth",1.25,"Color",'r',LineStyle='-');
 ylabel("log|\Psi (t,r*)|");
%plot(tarr_1(4:end),cell2mat(psi_extr_1(4:end)),"LineWidth",2,"Color",'r',LineStyle='-');
%ylabel("\Psi (t,r*)");
 xlabel("t")
 t=sprintf("N=%d, subd=%d, L=%d, FinalTime=%d, signal extracted at rstar=%1.4f",N,subd,L,FinalTime,x(r(1),c(1)));
 title(t)
 hold on

[psiarr_2,piarr_2,phiarr_2,tarr_2] = ScalarWaveKerr1D(psi_in,pi_in,phi_in,L,M,a,e(2),FinalTime);
tarr_2=cell2mat(tarr_2);

psi_extr_2={};
extr_posn = 250;
[r,c]=find(abs(x-extr_posn)<=1);
for i=1:length(tarr_2)
   psi_t=cell2mat(psiarr_2(i));
   psi=psi_t(r(1),c(1));
   psi_extr_2=[psi_extr_2,psi]; %#ok<AGROW> 
end
 psi_extr_2= cell2mat(psi_extr_2);

 plot(tarr_2(4:end),log10(abs(psi_extr_2(4:end))),"LineWidth",1.25,"Color",'b','LineStyle','--');
 ylabel("log|\Psi (t,r*)|");
%plot(tarr_2(4:end), cell2mat(psi_extr_2(4:end)),"LineWidth",2,"Color",'b','LineStyle','--');
%ylabel("\Psi (t,r*)");
xlabel("t")
t=sprintf("N=%d, subd=%d, L=%d, FinalTime=%d, signal extracted at rstar=%1.4f",N,subd,L,FinalTime,x(r(1),c(1)));
title(t)

hold off