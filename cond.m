function q = cond(par, k, T, A, idx, d, dir, flag)
    % Extract indices for z, r, and theta
    z_idx = idx{1}; 
    r_idx = idx{2}; 
    th_idx = idx{3};
    
    % Pre-allocate flux with the same size as the input indices
    q = zeros(size(z_idx));
    
    % Switch based on the direction flag
    switch dir
        case 1  % Axial direction (dz)
            % Axial flux calculation using d = dz
            if flag == 0
                % Flux in (to the LHS)
                q = k*A(z_idx, r_idx) .* ...
                    (T(z_idx-1, r_idx, th_idx) - T(z_idx, r_idx, th_idx)) ./ d;
            elseif flag == 1
                % Flux out (to the RHS)
                q = k*A(z_idx+1, r_idx) .* ...
                    (T(z_idx+1, r_idx, th_idx) - T(z_idx, r_idx, th_idx)) ./ d;
            end
            
        case 2  % Radial direction (dr)
            % Radial flux calculation using d = dr
            if flag == 0
                % Flux in (to the LHS)
                q = k * A(z_idx, r_idx) .* ...
                    (T(z_idx, r_idx-1, th_idx) - T(z_idx, r_idx, th_idx)) ./ d;
            elseif flag == 1
                % Flux out (to the RHS)
                q = k * A(z_idx, r_idx+1) .* ...
                    (T(z_idx, r_idx+1, th_idx) - T(z_idx, r_idx, th_idx)) ./ d;
            end

        case 3  % Angular direction (rdth)
            % Angular flux calculation using d = rdth(z, r)
            if flag == 0
                % Flux in (to the LHS)
                prev_th = mod(th_idx - 1 - 1, par.N_th) + 1; % ensure periodicity (modulo indexing)
                q = k*A(z_idx, r_idx) .* ...
                    (T(z_idx, r_idx, prev_th) - T(z_idx, r_idx, th_idx)) ./ d(z_idx, r_idx);
            elseif flag == 1
                % Flux out (to the RHS)
                next_th = mod(th_idx + 1 - 1, par.N_th) + 1; % ensure periodicity (modulo indexing)
                q = k*A(z_idx, r_idx) .* ...
                    (T(z_idx, r_idx, next_th) - T(z_idx, r_idx, th_idx)) ./ d(z_idx, r_idx);
            end
    end
end