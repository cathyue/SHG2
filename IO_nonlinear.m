ko1 = 3e6; %Hz
ke1 = 3e6; %Hz
ko2 = 6e6; %Hz
ke2 = 6e6; %Hz
delt = 1e9; %Hz
dw = 0;

c = 299792458; %m/s
R= 57e-6;
s0 = sqrt(1e-3);    %input 1mW
abs2g = 2.85e12*c/(2*pi*R);    %g(kmm), related to kappa_cal.m

fun = @NL_eq;
theta = pi/4;
x0 = sqrt(c/(2*pi*R)*ke1/((ko1+ke1)/2)^2*s0).*[cos(theta), sin(theta)];

a1 = fsolve(fun, x0);
abs2aout = ke2*abs2g*(a1(1)^2+a1(2)^2)^2/((delt)^2+(ko2+ke2)^2/4)
