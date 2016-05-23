function [ z ] = hl( l,x )
%UNTITLED Summary of this function goes here
%   spherical Hankel function

z = sqrt(pi./(2.*x)).*besselh(l+0.5, x);

end

