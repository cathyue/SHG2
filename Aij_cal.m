function [ Aij ] = Aij_cal( R, k1, n1, l1, k2, n2, l2 )
%UNTITLED2 此处显示有关此函数的摘要
%   my notebook, calculate A11, A12...

epsi0 = 8.854188e-12;   %F/m
Int1_1 = Int1_cal(R, k1, n1, l1);
Int1_2 = Int1_cal(R, k2, n2, l2);

funin = @(r) ((abs(phil(l1, n1.*k1.*r).*l1*(l1+1)./r.^2).^2....
+abs(l1.*n1.*k1./r.*(phil(l1-1, n1.*k1.*r)-l1./(n1.*k1.*r).*phil(l1, n1.*k1.*r))).^2)...
.*(abs(phil(l2, n2.*k2.*r).*l2*(l2+1)./r.^2).^2....
+abs(l2.*n2.*k2./r.*(phil(l2-1, n2.*k2.*r)-l2./(n2.*k2.*r).*phil(l2, n2.*k2.*r))).^2)).*r.^2;

A1 = n1^2*phil(l1, n1.*k1.*R)./kail(l1, k1.*R);
A2 = n2^2*phil(l2, n2.*k2.*R)./kail(l2, k2.*R);

funout = @(r) (abs(A1.*kail(l1, k1.*r).*l1*(l1+1)./r.^2).^2....
+abs(l1.*k1./r.*(A1.*kail(l1-1, k1.*r)-l1./(k1.*r).*A1.*kail(l1,k1.*r))).^2)...
.*(abs(A2.*kail(l2, k2.*r).*l2*(l2+1)./r.^2).^2....
+abs(l2.*k2./r.*(A2.*kail(l2-1, k2.*r)-l2./(k2.*r).*A2.*kail(l2, k2.*r))).^2).*r.^2;

IntS = integral(funin, 0, R) + integral(funout, R, R+2*pi/k1*20);

Aij = sqrt(l1.*l2).*2.*IntS.*sqrt(pi./((l1+l2)))./(pi^2.*epsi0^2.*n1.^2.*n2.^2.*Int1_1.*Int1_2);
end

