% IO_K_Thm.m
% This is to calculate the nonlinear coupled mode equations involving Kerr
% Thermal nonlinearities
para;

% tunable:
dw = 0;
s_in = sqrt(8.3e-4);    %input sqrt(W)
ke = [3e6; 6e6].*(2*pi);    %rad/s
ke0 = ke;
ke_s = 1;

% material para:
ko = [3e6; 6e6].*(2*pi);%rad/s, corresponding to Q~2e8
Q = 2*pi*c0./(lam0.*(ko+ke)./(2*pi));

P2m = zeros(1, length(ke_s));
P1m = zeros(1, length(ke_s));

for kke = 1:length(ke_s)
    ke = ke0.*ke_s(kke);
    for ks0 = 1:length(s_in)
        s0 = s_in(ks0);
        a0 = sqrt(ke).*s0./((ko+ke)./2);
        a0 = a0.*[1.1; 0];
        a0(2) = kmm*a0(1)^2/((ko(2)+ke(2))/2);
        Pext = [8.2132e-04, 4.5e-6];
%         kmm = 0;
        sweep = 9.1e10:0.5e8:9.5e10;
%         B = 0.*B;
%         delt = 0;
        theta0 = [0; pi/2];
        a0 = a0.*[1*(cos(theta0(1))+1i*sin(theta0(1))); 0.5*(cos(theta0(2))+1i*sin(theta0(2)))];
        x = [real(a0(1)), imag(a0(1)), real(a0(2)), imag(a0(2))];
        P1 = zeros(length(sweep),1);
        P2 = zeros(length(sweep),1);
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
            elseif (abs(func_final(3)^2+func_final(4)^2*ke(2))>=Pext(2))
                x0 = [1,1,0,0].*x;
                [func_final, fval, exitflag, output] = fsolve(@(x) NL_K_Thm(x, kmm, delt, dw, ke, ko, B, s0), x0, options);
            end
            flag(kw) = exitflag;
            a = [func_final(1)+1i*func_final(2); func_final(3)+1i*func_final(4)];
            Pout = abs(a.*sqrt(ke)).^2;
            P1(kw) = Pout(1);
            P2(kw) = Pout(2);
        end
        P2(P1>=Pext(1)) = 0;
        P1(P1>=Pext(1)) = 0;
        figure; plot(sweep, P1*1e-3, sweep, P2);
        ylim([0, 50e-7]);xlim([9e10, 9.4e10]);
    end
    [P2m(kke), ind] = max(P2);
        P1m(kke) = P1(ind);
        
end
% figure; plot(sweep, P1);