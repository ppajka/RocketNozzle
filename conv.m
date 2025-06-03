function q = conv(h, T, A, idx, T_inf, dir, flag)
    % Extract indices for z, r, and theta
    z_idx = idx{1}; 
    r_idx = idx{2}; 
    th_idx = idx{3};

    % Pre-allocate flux with the same size as the input indices
    q = zeros(size(T(z_idx, r_idx, th_idx)));

    % Switch based on the direction flag
    switch dir
        case 1  % Axial direction (dz)
            if flag == 0
                % Flux in (to the LHS)
                q = h .* A(z_idx, r_idx) .* ...
                    (T_inf - T(z_idx, r_idx, th_idx));
            else
                % Flux out (to the RHS)
                q = h .* A(z_idx+1, r_idx) .* ...
                    (T_inf - T(z_idx, r_idx, th_idx));
            end

        case 2  % Radial direction (dr)
            if flag == 0
                % Flux in (to the LHS)
                q = h .* A(z_idx, r_idx) .* ...
                    (T_inf - T(z_idx, r_idx, th_idx));
            else
                % Flux out (to the RHS)
                q = h .* A(z_idx, r_idx+1) .* ...
                    (T_inf - T(z_idx, r_idx, th_idx));
            end
    end
end