% Calculate the bulk quadruaple term delta(E.*D)E
clear; clc;
para;

epsi0 = 8.854188e-12;   %F/m
deltaP = 7.8e-22;   %m2/V

Int1_1 = Int1_cal(R, k0(1), n0(1), l0(1));
Int1_2 = Int1_cal(R, k0(2), n0(2), l0(2));

funr_in = @(r) conj(phil(l0(2), n0(2).*k0(2).*r)).*phil(l0(1), n0(1).*k0(1).*r)./r.^4.*...
    [phil(l0(1), n0(1).*k0(1).*r).*(-2./r-l0(1)./(r))+phil(l0(1)-1, n0(1).*k0(1).*r).*n0(1).*k0(1)];
A1 = n0(1)^2*phil(l0(1), n0(1).*k0(1).*R)./kail(l0(1), k0(1).*R);
A2 = n0(2)^2*phil(l0(2), n0(2).*k0(2).*R)./kail(l0(2), k0(2).*R);
funr_out = @(r) conj(A2.*kail(l0(2), k0(2).*r)).*A1.*kail(l0(1), k0(1).*r)./r.^4.*...
    [A1.*kail(l0(1), k0(1).*r).*(-2./r-l0(1)./(r))+A1.*kail(l0(1)-1, k0(1).*r).*k0(1)];

Int_r = integral(funr_in, 0, R)+integral(funr_out, R, R+2*pi/k0(1)*20);

Coef = 2.*1i.*w0(1).^2./w0(2)./n0(2).^2.*deltaP.*...
    2^0.25.*l0(1).^(3/4)./(pi.^(9/4).*epsi0.^(3/2).*n0(1).^2.*n0(2).*Int1_1.*sqrt(Int1_2)).*2.*l0(1).^3.*(2.*l0(1)+1).*(l0(1)+1).^2.*...
    sqrt(pi./2./l0(1)).*2.*pi./(2./epsi0./n0(2).^2).*Int_r;
% da2/dt=Coef*a1^2; checked170124

%% for  comparison, the surface term

Coef_s = 2.*1i.*w0(1).^2./w0(2)./n0(2).^2.*R.^2.*kai_ttt.*phil(l0(2),n0(2).*k0(2).*R).*2*l0(1).*(2*l0(1)+1).*...
    phil(l0(1), n0(1).*k0(1).*R).^2.*l0(1).^2.*(l0(1)+1).^2./R^6.*2^0.25.*l0(1)^0.75./...
    (pi^(9/4).*epsi0^(3/2).*n0(1)^2.*n0(2).*sqrt(Int1_2).*Int1_1).*sqrt(pi/2/l0(1)).*2.*pi./(2./epsi0./n0(2).^2);

%% TE+TE=TM, bulk term, but use the eigen frequencies from TM

S_in1 = @(r) abs(phil(l0(1), n0(1).*k0(1).*r)).^2;
A1te = phil(l0(1), n0(1).*k0(1).*R)./kail(l0(1), k0(1).*R);
S_out1 = @(r) abs(A1te.*kail(l0(1), k0(1).*r)).^2;
IntS1 = integral(S_in1, 0, R)+integral(S_out1, R, R+2*pi/k0(1)*20);

S_in2 = @(r) abs(phil(2.*l0(1)-1, n0(2).*k0(2).*r)).^2;
A2te = phil(2.*l0(1)-1, n0(2).*k0(2).*R)./kail(2.*l0(1)-1, k0(2).*R);
S_out2 = @(r) abs(A2te.*kail(2*l0(1)-1, k0(2).*r)).^2;
IntS2 = integral(S_in2, 0, R)+integral(S_out2, R, R+2*pi/k0(1)*20);

funte_in = @(r) conj(phil(2.*l0(1)-1, n0(2).*k0(2).*r)).*phil(l0(1), n0(1).*k0(1).*r).^2./r.^2;
funte_out = @(r) conj(A2te.*kail(2*l0(1)-1, k0(2).*r)).*(A1te.*kail(l0(1), k0(1).*r)).^2./r.^2;
Intte = integral(funte_in, 0, R)+integral(funte_out, R, R+2*pi/k0(1)*20);

Coef_te = 2.*1i.*w0(1).^2./w0(2)./n0(2).^2.*deltaP.*...
    2.*pi.*Intte.*1i./(epsi0.^1.5.*n0(1).^2.*n0(2).*(2*pi).^1.5.*IntS1.*sqrt(IntS2)).*(l0(1)-1)./(4.*l0(1)-3).*...
    (2.*l0(1)-2).^1.5.*(2./pi).^0.25.*(4.*l0(1)-4).^0.75./sqrt(4.*l0(1)-6)./(2./epsi0./n0(2).^2);

%% double check
C_TM = Int_r.*l0(1).^0.25./(2.^0.25);

fte_in = @(r) conj(phil(2.*l0(1)-1, n0(2).*k0(2).*r)).*phil(l0(1), n0(1).*k0(1).*r).^2./r.^5;
fte_out = @(r) conj(A2te.*kail(2*l0(1)-1, k0(2).*r)).*(A1te.*kail(l0(1), k0(1).*r)).^2./r.^5;
Int = integral(fte_in, 0, R)+integral(fte_out, R, R+2*pi/k0(1)*20);
C_TE = Int.*l0(1).^0.75./(8.^0.75);