function dT_dt_flu = flu_mean_temp(z, noz, idx, flu, chnl, par)
    % Extract index for z
    z_idx = idx; 

    % Extract parameters
    A_r = par.A_r;
    A_th = par.A_th;
    dz = par.dz(2);

    % Extract properties
    rho = flu.rho;
    m_dot = flu.m_dot;
    cp = flu.cp;

    % Wetted area at perimeter of channel
    A_conv = (A_r(:, chnl.r_i_idx)*(chnl.th_o_idx - chnl.th_i_idx + 1) ...  % radial areas at chnl.r_i
        + A_r(:, chnl.r_o_idx+1)*(chnl.th_o_idx - chnl.th_i_idx + 1) ...  % radial areas at chnl.r_o
        + 2*sum(A_th(:, chnl.r_i_idx:chnl.r_o_idx), 2));  % circumferential areas at both chnl.th_i and chnl.th_o

    % Use upwind:
    if z_idx == 1
        T_up = flu.T_in;
    else
        T_up = flu.T(z_idx-1);
    end

    % Combine advection and convection terms to find dT/dt of the fluid
    A_flu = chnl.A_flu_func(z(z_idx));
    A_flu = A_flu(:);  % make sure this area is a col vec
    dT_dt_adv = m_dot ./ (rho*A_flu*dz) ...
        .* (T_up - flu.T(z_idx));
    dT_dt_conv = flu.h_T(z_idx) .* A_conv(z_idx) ./ (rho*cp*A_flu*dz) ...
        .* (noz.T_m(z_idx) - flu.T(z_idx));
    dT_dt_flu = dT_dt_adv + dT_dt_conv;

end