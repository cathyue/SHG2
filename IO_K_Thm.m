% IO_K_Thm.m
% This is to calculate the nonlinear coupled mode equations involving Kerr
% Thermal nonlinearities
para;   %warning: You should change sample if you use different pairs of paras!!!

% tunable:
dw = 0;
theta_sin = 0;
s_in = sqrt([3.60:0.1:4.9].*1e-4).*(cos(theta_sin)+1i*sin(theta_sin));    %input sqrt(W)
ke = [3e6; 6e6].*(2*pi);    %rad/s
ke0 = ke;
ke_s = 1.6; %coupling coeffecient, ke/ko

% material para:
ko = [3e6; 6e6].*(2*pi);%rad/s, corresponding to Q~2e8
Q = 2*pi*c0./(lam0.*(ko+ke)./(2*pi));

P2m = zeros(1, length(s_in));
P1m = zeros(1, length(s_in));
P2id = zeros(1, length(s_in));
P1id = zeros(1, length(s_in));
B0 = B;
delt0 = delt;

for kke = 1:length(ke_s)
    ke = ke0.*ke_s(kke);
    for ks0 = 1:length(s_in)
        s0 = s_in(ks0);
        a0 = sqrt(ke).*abs(s0)./((ko+ke)./2);
        a0 = a0.*[1.1; 0];
        a0(2) = kmm*a0(1)^2/((ko(2)+ke(2))/2);
%         Pext = [4.9686e-04, 15.678e-7];
                Pext = [1,1];
        %         kmm = 0;
        sweep = 3.75e10:0.01e8:3.9e10;%-1e8:1e6:1e8;%
        B = B0;
        delt = delt0;
        theta0 = [0; pi/2];
        a0 = a0.*[1*(cos(theta0(1))+1i*sin(theta0(1))); 0.5*(cos(theta0(2))+1i*sin(theta0(2)))];
        x = [real(a0(1)), imag(a0(1)), real(a0(2)), imag(a0(2))];
        P1 = zeros(length(sweep),1);
        P2 = zeros(length(sweep),1);
        trans = zeros(length(sweep),1);
        flag = zeros(length(sweep), 1);
        options = optimoptions('fsolve','MaxFunEval', 400,'MaxIter', 400, 'algorithm', 'levenberg-marquardt');
        for kw = 1:length(sweep)
            dw = sweep(kw);
            %     if kw<84
            %         x0 = 0.*x;
            %     else
            %         x0 = x;
            %     end
            x0 = x;
            %     exitflag = -1;
            [func_final, fval, exitflag, output] = fsolve(@(x) NL_K_Thm(x, kmm, delt, dw, ke, ko, B, s0), x0, options);
            
            if (abs(func_final(1)^2+func_final(2)^2)*ke(1)>=Pext(1))...
                    ||(exitflag<0)
                x0 = 0.*x;
                [func_final, fval, exitflag, output] = fsolve(@(x) NL_K_Thm(x, kmm, delt, dw, ke, ko, B, s0), x0, options);
            elseif ((abs(func_final(3)^2+func_final(4)^2)*ke(2))>=Pext(2))
                x0 = [1,1,0,0].*x;
                [func_final, fval, exitflag, output] = fsolve(@(x) NL_K_Thm(x, kmm, delt, dw, ke, ko, B, s0), x0, options);
            end
            flag(kw) = exitflag;
            a = [func_final(1)+1i*func_final(2); func_final(3)+1i*func_final(4)];
            Pout = abs(a.*sqrt(ke)).^2;
            P1(kw) = Pout(1);
            P2(kw) = Pout(2);
            trans(kw) = -s0+sqrt(ke(1))*a(1);
        end
        P2(P1>=Pext(1)) = 0;
        P1(P1>=Pext(1)) = 0;
        figure; plot(sweep, P1*1e-3, sweep, P2);
        ylim([0, 8e-7]);xlim([sweep(1), sweep(length(sweep))]);
        [P2m(ks0), ind] = max(P2);
        P1m(ks0) = P1(ind);
        dw = 0;
        B = 0.*B0;
        delt = 0;
        x0 = x;
        [func_final, fval, exitflag, output] = fsolve(@(x) NL_K_Thm(x, kmm, delt, dw, ke, ko, B, s0), x0, options);
        a = [func_final(1)+1i*func_final(2); func_final(3)+1i*func_final(4)];
        Pout = abs(a.*sqrt(ke)).^2;
        P1id(ks0) = Pout(1);
        P2id(ks0) = Pout(2);
    end
    
end
% figure; plot(sweep, P1);