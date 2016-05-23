function [ F ] = NL_eq( a1 )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明

c = 299792458; %m/s
R = 57e-6;  %m/s

ko1 = 3e6; %Hz
ke1 = 3e6; %Hz
ko2 = 6e6; %Hz
ke2 = 6e6; %Hz
delt = 1e9; %Hz
kap1 = (ko1+ke1)/2;
kap2 = (ko2+ke2)/2;
dw = 0;
abs2g = 2.85e12*c/(2*pi*R);    %g(kmm), related to kappa_cal.m
coef_g = abs2g/((2*dw+delt)^2+(kap2)^2);

s = sqrt(1e-3); %sqrt(Watt)

% F = -a1*(1i*dw+(ko1+ke1)/2)-abs2g*a1*abs(a1)^2/(1i*(2*dw+delt)+(ko2+ke2)/2)+sqrt(ke1)*s;

F(1) = -(a1(1)*kap1-a1(2)*dw)-coef_g*(a1(1)^2+a1(2)^2)*(a1(1)*kap2+a1(2)*(2*dw+delt))+sqrt(ke1)*s;
F(2) = -(dw*a1(1)+a1(2)*kap1)-coef_g*(a1(1)^2+a1(2)^2)*(a1(2)*kap2-a1(1)*(2*dw+delt));


end

