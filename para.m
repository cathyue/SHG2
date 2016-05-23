%   para.m
load sample3.mat;

% from sample.mat, close to on resonance condition
lam0 = [lam10; lam20];
n0 = [n10; n20];
w0 = [w10; w20];
l0 = [l0; 2*l0];
delt = w0(2)-w0(1)*2; %from detuning.m

% constant:
c0 = 299792458; %m/s
Z0 = 376.73;    %free space resistance
epsi0 = 8.854188e-12;   %F/m

%Kerr nonlinearity
n2 = [2.79e-20; 2.48e-20]; %m2/W, Review and assessment of measured values of the nonlinear refractive-index coefficient of fused silica David Milam
kai3 = n2.*4.*n0.^2.*epsi0*c0/3;   %m2/V2

% Thermal nonlinearity
dndT = 6e-6;  %1/K
rho = 2200; %kg/m3
C = 740;    %J/(kgK)
D = 9.5e-7; %m2/s
Qab = [7e8; 2e10];  %from Rokhsari_APL_2004, Fig.2
b = [1.7e-6;sqrt(1.08^2+1)*1e-6];    %m, estimated from phil(l,x)
dthet = 2*D./(b.^2)./124.54/2;    %emperical value

%Bij_cal
B = zeros(2,2);
for i = 1:2
    for j = 1:2
        Aij = Aij_cal(R, 2*pi/(lam0(i)), n0(i), l0(i), 2*pi/(lam0(j)), n0(j), l0(j));
        
        same_ = (i==j);
        
        B(i, j) = (3*(1+same_)*kai3(j)*w0(j)*Aij/n0(j)^2 ...
            + epsi0*w0(j)/n0(j)*dndT/(rho*C*dthet(j))*n0(i)^2*w0(i)/Qab(i)*Aij)/(2/(epsi0*n0(j)^2));
    end
end

% kappa & g cal
kai_ttt = 59e-22; % m2/V, second order susceptibility, surface effective
kai_tll = 3.8e-22;
kai_llt = 7.9e-22;
k0 = 2*pi./lam0;
zl = (1+n0(1))./2.*jl(l0(1), n0(1).*k0(1)*R);
drzl_dr = 0.5.*(n0(1).*k0(1).*R.*jl(l0(1)-1, n0(1).*k0(1).*R) ...
            -l0(1).*(1+n0(1).^2).*jl(l0(1),n0(1).*k0(1).*R) ...
            +n0(1).^2.*k0(1)*R.*jl(l0(1), n0(1).*k0(1).*R).*hl(l0(1)-1, k0(1).*R)./hl(l0(1), k0(1).*R));
Gmm = Gm2(l0(2),k0(2)*R,n0(2));

kmm = c0*n0(1)^2*k0(1)^4*zl^2*l0(1)/(sqrt(2*epsi0)*n0(2)*R*Gmm*k0(2)) ...
    *sqrt(Int1_cal(R, k0(2), n0(2), l0(2)))/Int1_cal(R, k0(1), n0(1), l0(1)) ...
    *(kai_ttt-1/(l0(1)*zl)^2*(drzl_dr)^2*kai_tll) ...
    *(-1)^(2*l0(1)+l0(2))*sqrt(l0(1))*l0(2)^0.25/(sqrt(2)*pi^0.75*sqrt(l0(1)+0.5*l0(2)));

% % sq cal, to change from kmm (PRA2008) to g
% sq = [0;0];
% for ksq = 1:2
%     sq(ksq) = k0(ksq)*l0(ksq)*phil(l0(ksq), n0(ksq).*k0(ksq).*R)*sqrt(c0) ...
%         /(sqrt(Int1_cal(R, k0(ksq), n0(ksq), l0(ksq)))*kail(l0(ksq), k0(ksq).*R));
% end
% 
% g = [0;0];
% g(1) = conj(kmm)*sq(2)*conj(sq(1))/sq(1);
% g(2) = kmm*sq(1)^2/sq(2);



