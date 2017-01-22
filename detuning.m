 % detuning between the 1 and 2 order wave, regard to different l (x-axis),
 % and p (different lines)
l = 60:200;
p = 1;
n0 = 1.44;  %initial refractive index 
N0 = 1.45;
P = [1,2,3];
R = 60e-6/2;
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
coeP2 = zeros(length(P), length(l));
c0 = 299792458; %m/s
lam1 = (2*pi*c0./w1).';
lam2 = (2*pi*c0./w2).';
Q = 2e8;

figure;hold on;
for kp = 1:length(P)
    dw(kp, :) = (w2(kp, :)-2.*w1)./w2(kp, :);
    coeP2(kp,:) = (1./(1+(2.*dw(kp,:).*Q).^2)).';
    plot(lam1, dw(kp, :));
end
% semilogy(lam1, coeP2(:),'o');
coeP2 = coeP2.';

%legend('P=1', '2','3','4');
xlabel('lam1');ylabel('CoeP2');
plot (lam1, zeros(1,length(l)));


   