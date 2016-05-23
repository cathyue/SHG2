function [ z ] = phil( l,x )
%UNTITLED2 Summary of this function goes here
%   spherical Bessel function

z = sqrt(x.*pi./(2)).*besselj(l+0.5, x);


end

