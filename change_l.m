kl0 = 169:174;
lam1s = zeros(length(kl0),1);
delts = zeros(length(kl0),1);
Pins = zeros(length(kl0),1);
P2s = zeros(length(kl0),1);
for kkl0 = 1:length(kl0)
    w10 = w1(find(l==kl0(kkl0)));
w20 = w2(find(l==kl0(kkl0)));
c0 = 299792458; %m/s
lam10 = 2*pi*c0/w10;
lam1s(kkl0) = lam10;
lam20 = 2*pi*c0/w20;
n10 = n_lam(lam10*1e6);
n20 = n_lam(lam20*1e6);
% save('samp31_171', 'lam10','lam20','w10','w20','n10','n20','l0','R');

% from sample.mat, close to on resonance condition
lam0 = [lam10; lam20];
n0 = [n10; n20];
w0 = [w10; w20];
l0 = [kl0(kkl0); 2*kl0(kkl0)];
delt = w0(2)-w0(1)*2; %from detuning.m
delts(kkl0) = delt;

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




% tunable:
dw = 0;
ke = [3e6; 6e6].*(2*pi);    %rad/s
ke0 = ke;
ke_s = 1; %coupling coeffecient, ke/ko

% material para:
ko = [3e6; 6e6].*(2*pi);%rad/s, corresponding to Q~2e8
Q = 2*pi*c0./(lam0.*(ko+ke)./(2*pi));

%estimate the power to ensure SHG
a12_est = delt/(B(1,2)-2*B(1,1));
wc_est = a12_est*B(1,1);
Pin_est = a12_est*(ke(1)+ko(1))^2/4/ke(1);
Pins(kkl0) = Pin_est+0.1e-4;
theta_sin = 0;
s_in = sqrt(Pin_est+0.1e-4).*(cos(theta_sin)+1i*sin(theta_sin));    %input sqrt(W)

P2m = zeros(1, length(s_in));
P1m = zeros(1, length(s_in));
P2id = zeros(1, length(s_in));
P1id = zeros(1, length(s_in));
B0 = B;
delt0 = delt;

for kke = 1:length(ke_s)
    ke = ke0.*ke_s(kke);
    for ks0 = 1:length(s_in)
        s0 = s_in(ks0);
        a0 = sqrt(ke).*abs(s0)./((ko+ke)./2);
        a0 = a0.*[1.1; 0];
        a0(2) = kmm*a0(1)^2/((ko(2)+ke(2))/2);
%         Pext = [3.3005e-04, 9.0811e-7];  % change manually!
%         Pext = [1,1];
        theta0 = [0; pi/2];
        a0 = a0.*[1*(cos(theta0(1))+1i*sin(theta0(1))); 0.5*(cos(theta0(2))+1i*sin(theta0(2)))];
        x = [real(a0(1)), imag(a0(1)), real(a0(2)), imag(a0(2))];
        options = optimoptions('fsolve','MaxFunEval', 400,'MaxIter', 400, 'algorithm', 'levenberg-marquardt');
        
        
%         sweep = wc_est-0.3e10:0.05e8:wc_est+0.3e10;%9.2e10:0.01e8:9.4e10;%-1e8:1e6:1e8;%
%         B = B0;
%         delt = delt0;
%         P1 = zeros(length(sweep),1);
%         P2 = zeros(length(sweep),1);
%         trans = zeros(length(sweep),1);
%         flag = zeros(length(sweep), 1);
%         for kw = 1:length(sweep)
%             dw = sweep(kw);
%             %     if kw<84
%             %         x0 = 0.*x;
%             %     else
%             %         x0 = x;
%             %     end
%             x0 = x;
%             %     exitflag = -1;
%             [func_final, fval, exitflag, output] = fsolve(@(x) NL_K_Thm(x, kmm, delt, dw, ke, ko, B, s0), x0, options);
%             
%             if (abs(func_final(1)^2+func_final(2)^2)*ke(1)>=Pext(1))...
%                     ||(exitflag<0)
%                 x0 = 0.*x;
%                 [func_final, fval, exitflag, output] = fsolve(@(x) NL_K_Thm(x, kmm, delt, dw, ke, ko, B, s0), x0, options);
%             elseif ((abs(func_final(3)^2+func_final(4)^2)*ke(2))>=Pext(2))
%                 x0 = [1,1,0,0].*x;
%                 [func_final, fval, exitflag, output] = fsolve(@(x) NL_K_Thm(x, kmm, delt, dw, ke, ko, B, s0), x0, options);
%             end
%             flag(kw) = exitflag;
%             a = [func_final(1)+1i*func_final(2); func_final(3)+1i*func_final(4)];
%             Pout = abs(a.*sqrt(ke)).^2;
%             P1(kw) = Pout(1);
%             P2(kw) = Pout(2);
%             trans(kw) = -s0+sqrt(ke(1))*a(1);
%         end
%         P2(P1>=Pext(1)) = 0;
%         P2(P2>=Pext(2)) = 0;
%         P1(P1>=Pext(1)) = 0;
%         figure; plot(sweep, P1*1e-3, sweep, P2);
%         ylim([0, Pext(2)]);xlim([sweep(1), sweep(length(sweep))]);
%         [P2m(ks0), ind] = max(P2);
%         P1m(ks0) = P1(ind);
        
                dw = 0;
                B = 0.*B0;
                delt = 0;
                x0 = x;
                [func_final, fval, exitflag, output] = fsolve(@(x) NL_K_Thm(x, kmm, delt, dw, ke, ko, B, s0), x0, options);
                a = [func_final(1)+1i*func_final(2); func_final(3)+1i*func_final(4)];
                Pout = abs(a.*sqrt(ke)).^2;
                P1id(ks0) = Pout(1);
                P2id(ks0) = Pout(2);
                P2s(kkl0) = Pout(2);
    end
    
end
end