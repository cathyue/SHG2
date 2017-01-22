clear; clc;

%% Basic parameters
c0 = 299792458;

%% Fundamental wave
L0 = (162:168)+5;
p = 1:10;
n0 = 1.44;  %initial refractive index
N0 = 1.45;
% P = 2;
R = 61e-6/2;
w1TM = zeros(length(p), length(L0));
w1TE = zeros(length(p), length(L0));

l = [L0;(L0-8);(L0-14);(L0-20);(L0-24);...
    (L0-28);(L0-33);(L0-37);(L0-42);(L0-46)];

%TM
Plr = 0;
for kl = 1:length(L0)
    for kp = 1:length(p)
%         Lin = l(kl)-(p(kp)-1)*6;
        [w1TM(kp,kl), n] = ome_lpP(l(kp,kl),p(kp),Plr, n0, R);
        [w1TM(kp,kl), n] = ome_lpP(l(kp,kl),p(kp),Plr, n, R);
        [w1TM(kp,kl), n] = ome_lpP(l(kp,kl),p(kp),Plr, n, R);
        [w1TM(kp,kl), n] = ome_lpP(l(kp,kl),p(kp),Plr, n, R);
    end
end

%TE
Plr = 1;
for kl = 1:length(L0)
    for kp = 1:length(p)
        [w1TE(kp,kl), n] = ome_lpP(l(kp,kl),p(kp),Plr, n0, R);
        [w1TE(kp,kl), n] = ome_lpP(l(kp,kl),p(kp),Plr, n, R);
        [w1TE(kp,kl), n] = ome_lpP(l(kp,kl),p(kp),Plr, n, R);
        [w1TE(kp,kl), n] = ome_lpP(l(kp,kl),p(kp),Plr, n, R);
    end
end

lam1TM = c0*2*pi./w1TM.*1e9;    %unit: nm
lam1TE = c0*2*pi./w1TE.*1e9;    %unit: nm

%% Second Harmonic
l = 326:338;
p = 1:10;
n0 = 1.45;
% P = 2;
R = 61e-6/2;
w2TM = zeros(length(p), length(l));
w2TE = zeros(length(p), length(l));

%TM
Plr = 0;
for kl = 1:length(l)
    for kp = 1:length(p)
        [w2TM(kp,kl), n] = ome_lpP(l(kl),p(kp),Plr, n0, R);
        [w2TM(kp,kl), n] = ome_lpP(l(kl),p(kp),Plr, n, R);
        [w2TM(kp,kl), n] = ome_lpP(l(kl),p(kp),Plr, n, R);
        [w2TM(kp,kl), n] = ome_lpP(l(kl),p(kp),Plr, n, R);
    end
end

%TE
Plr = 1;
for kl = 1:length(l)
    for kp = 1:length(p)
        [w2TE(kp,kl), n] = ome_lpP(l(kl),p(kp),Plr, n0, R);
        [w2TE(kp,kl), n] = ome_lpP(l(kl),p(kp),Plr, n, R);
        [w2TE(kp,kl), n] = ome_lpP(l(kl),p(kp),Plr, n, R);
        [w2TE(kp,kl), n] = ome_lpP(l(kl),p(kp),Plr, n, R);
    end
end

lam2TM = c0*2*pi./w2TM.*1e9;    %unit: nm
lam2TE = c0*2*pi./w2TE.*1e9;    %unit: nm