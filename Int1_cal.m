function [ Int1 ] = Int1_cal( R, k, n, l )
%UNTITLED 此处显示有关此函数的摘要
%   my notebook, ->only taking into account m=l>>1

funin = @(r) l.^2.*(phil(l-1, n.*k.*r)-l./(n.*k.*r).*phil(l, n.*k.*r)).^2.*(n.*k).^2 ...
    +(phil(l, n.*k.*r).*l.*(l+1)./r).^2;
A = n^2*phil(l, n.*k.*R)./kail(l, k.*R);
funout = @(r) l.^2.*abs(A.*kail(l-1, k.*r)-l./(k.*r).*A.*kail(l, k.*r)).^2.*k.^2....
    +abs(A.*kail(l, k.*r).*l.*(l+1)./r).^2;

Int1 = integral(funin, 0, R)+integral(funout, R, R+2*pi/k*20);

%checked, 20160518

end

