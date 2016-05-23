function [ F ] = Fm( l, x, n )
%UNTITLED3 Summary of this function goes here
%   Kozyreff PRA 2008, eq. (34)

F = hl(l-1,x)-hl(l,x)*jl(l-1,x)/(n*jl(l,x))- (l+1)/(x)*hl(l,x)*(1-1/(n^2));

end

