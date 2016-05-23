function [ z ] = jl( l,x )
%UNTITLED2 Summary of this function goes here
%   spherical Bessel function

z = sqrt(pi./(2*x)).*besselj(l+0.5, x);


end

