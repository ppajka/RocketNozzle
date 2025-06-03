function [h_comb, T_comb] = Bartz(T_w, A)
    mu = 1.059E-4;
    cp = 3153.1;
    cv = 2673.8;
    gamma = cp/cv;
    mdot = 514.50;
    
    D_star = 0.539;
    A_star = (D_star/2)^2*pi;
    
    T_c = 3300+273.15;
    P_c = 20.643E5;
    
    c_star = P_c*A_star/mdot;
    
    Pr = 0.767;
    r_c = 0.5;

    A_r = A./A_star;
    %A_r(1) = 1.01;
    for i = 1:length(A)
        M(i) = Mach(A_r(i), gamma);
    end

    sigma = 1./((0.5.*T_w./T_c.*(1+(gamma-1)./2.*M.^2)+0.5).^(0.8-0.6./5).*(1+(gamma-1)./2.*M.^2).^0.6./5);

    h_comb = (0.026./D_star.^0.2).*((mu.^0.2.*cp)./(Pr.^0.6)).*((P_c./c_star).^0.8).*((D_star./r_c).^0.1).*((A_star./A).^0.9).*sigma;
    
    T_comb = (1+(gamma-1)./2.*M.^2).^-1.*T_c;
end

function M = Mach(A_r, g)
    exp = (g+1)/(2*(g-1));
    F = @(x) ((g+1)./2).^-exp.*(((1+(g-1)./2.*x.^2).^exp)./x) - A_r;
    options = optimset('Display','off');
    M = fsolve(F, 1.01, options);
end