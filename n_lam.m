function [ n ] = n_lam( lam )
% calculate n_lam(microns), Kozyreff Appendix D

B1 = 0.6961663;
B2 = 0.4079426;
B3 = 0.8974794;
L1 = 0.0684043; %microns
L2 = 0.1162414;
L3 = 9.896161;
n = sqrt(1+(lam^2)*(B1/(lam^2-L1^2)+B2/(lam^2-L2^2)+B3/(lam^2-L3^2)));

end

