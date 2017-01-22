function [ omega, n ] = ome_lpP( l,p,Plr,n, R )
%UNTITLED5 Summary of this function goes here
%   resonant frequence in a spherical resonator, p<=12, TM mode, n is the
%   initial approximate value
%   pg 9, equation 3, http://dx.doi.org/10.1364/AOP.7.000168   (sensor review, Vollmer)

airyroots = -[-2.338107410459764, -4.087949444130973, -5.520559828095556, -6.786708090071763, ...
    -7.944133587120849, -9.022650853340981, -10.040174341558084, -11.008524303733264,...
    -11.936015563236262, -12.828776752865757, -13.691489035210719, -14.527829951775335];
c = 299792458;
%n = n_lam(lam*1e6);
if Plr ==0  %TM
    P = 1/n;
else %TE
    P = n;
end
nu = l+0.5;

omega = c/(n*R)*(nu+airyroots(p)*nu^(1/3)/(2^(1/3))-P/sqrt(n^2-1)+0.3*airyroots(p)^2/(2^(2/3)*nu^(1/3))...
    -P*(n^2-2*P/3)*airyroots(p)/((n^2-1)^1.5*2^(1/3)*nu^(2/3)));
n1 = n_lam(2*pi*c/omega*1e6);
% n = n1;
% omega1 = c/(n*R)*(nu+airyroots(p)*nu^(1/3)/(2^(1/3))-P/sqrt(n^2-1)+0.3*airyroots(p)^2/(2^(2/3)*nu^(1/3))...
%     -P*(n^2-2*P/3)*airyroots(p)/((n^2-1)^1.5*2^(1/3)*nu^(2/3)));

while abs(n1-n) > 1e-15
    n = n1;
    omega = c/(n*R)*(nu+airyroots(p)*nu^(1/3)/(2^(1/3))-P/sqrt(n^2-1)+0.3*airyroots(p)^2/(2^(2/3)*nu^(1/3))...
        -P*(n^2-2*P/3)*airyroots(p)/((n^2-1)^1.5*2^(1/3)*nu^(2/3)));
    n1 = n_lam(2*pi*c/omega*1e6);
end





end

