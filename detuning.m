 % detuning between the 1 and 2 order wave, regard to different l (x-axis),
 % and p (different lines)
l = 60:300;
p = 1;
n0 = 1.44;  %initial refractive index 
N0 = 1.45;
P = 2;
R = 30e-6;
w1 = zeros(1, length(l));
w2 = zeros(length(P),length(l));

for kl = 1:length(l)
    [w1(kl), n] = ome_lp(l(kl),p, n0, R);
    [w1(kl), n] = ome_lp(l(kl),p, n, R);
    [w1(kl), n] = ome_lp(l(kl),p, n, R);
    [w1(kl), n] = ome_lp(l(kl),p, n, R);
    for kp = 1:length(P)
        [w2(kp, kl), n] = ome_lp(2*l(kl), P(kp), N0, R);
        [w2(kp, kl), n] = ome_lp(2*l(kl), P(kp), n, R);
        [w2(kp, kl), n] = ome_lp(2*l(kl), P(kp), n, R);
        [w2(kp, kl), n] = ome_lp(2*l(kl), P(kp), n, R);
    end
end

dw = zeros(length(P), length(l));

figure;
hold on;
for kp = 1:length(P)
    dw(kp, :) = (w2(kp, :)-2.*w1)./w2(kp, :);
    plot(l, dw(kp, :));
end

%legend('P=1', '2','3','4');
xlabel('l');ylabel('\Delta/\omega_2')
plot (l, zeros(1,length(l)));


   