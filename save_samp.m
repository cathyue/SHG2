l0 = 169;   %change!
w10 = w1(find(l==l0));
w20 = w2(find(l==l0));
c0 = 299792458; %m/s
lam10 = 2*pi*c0/w10;
lam20 = 2*pi*c0/w20;
n10 = n_lam(lam10*1e6);
n20 = n_lam(lam20*1e6);
save('samp30_169', 'lam10','lam20','w10','w20','n10','n20','l0','R');
%change