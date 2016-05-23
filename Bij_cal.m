function [ Bij ] = Bij_cal( i, j )
%UNTITLED 此处显示有关此函数的摘要
%   my notebook, Bij_cal

para

Aij = Aij_cal(R, 2*pi/(lam0(i)), n0(i), l0(i), 2*pi/(lam0(j)), n0(j), l0(j));

same_ = (i==j);

Bij = (3*(1+same_)*kai3(j)*w0(j)*Aij/n0(j)^2 ...
    +epsi0*w0(j)/n0(j)*dndT/(rho*C*dthet(j))*n0(i)^2*w0(i)/Qab(i)*Aij)/(2/(epsi0*n0(j)^2));


end

