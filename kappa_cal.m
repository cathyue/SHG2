% P2 = coeff*P1^2/(1+(dw*Q2/w2)); 
% This is to calculate coeff, unit: 1/Watt

load sample.mat
Z0 = 376.73;    %free space resistance
l = l0;
% p = 1;
L=2*l;
a= R;
lam1 = lam10; %m
lam2 = lam20;
n1 = n10;
n2 = n20;
epsi = 8.8541878176e-12; %F/m
kai_ttt = 59e-22; % m2/V, second order susceptibility, surface effective
kai_tll = 3.8e-22;
kai_llt = 7.9e-22;
k1 = 2*pi/lam1;
k2 = 2*pi/lam2;

% ko1 = 3e6; %Hz
% ke1 = 3e6; %Hz
% ko2 = 6e6; %Hz
% ke2 = 6e6; %Hz
% delt = 1e10; %Hz

% Q1 = 1e8;
% Q2 = 1e8;
c = 299792458; %m/s

zl = hl(l, k1*a);
% dr = a/(1e15);
drzl_dr = -l*hl(l,k1*a)+k1*a*hl(l-1, k1*a);
Gmm = Gm(L,k2*a,n2);

kmm = sqrt(2/Z0)*n1^2*k1^2*zl^2/(epsi*n2*a*Gmm)*(kai_ttt-1/(l*zl)^2*(drzl_dr)^2*kai_tll)...
    *sqrt(l/(4*pi))*0.5
% coeff = abs(kmm)^2*lam1^2*Q1^2*Q2*n2*lam2/(4*pi^4*a*n1^2*c^2)
% estimated coeff ~0.02
