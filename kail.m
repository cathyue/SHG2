function [ z ] = kail( l,x )
%UNTITLED Summary of this function goes here
%   spherical Hankel function

z = sqrt(x.*pi./(2)).*besselh(l+0.5, x);

end

