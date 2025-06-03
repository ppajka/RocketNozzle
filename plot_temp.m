function plot_temp(T_soln, noz, par, z, r, th, idx)
    % Extract the 2D slice
    T_slice = squeeze(T_soln(idx.z, idx.r, idx.th));
    
    % Determine which coordinate is fixed
    if ~isequal(idx.z, ':')
        is2D = 1;
        % Fixed z: (r, θ) plane
        [R, Th] = meshgrid(r, th);
        X = R;
        Y = Th;
        x_label = 'r [m]';
        y_label = 'θ [rad]';
        
    elseif ~isequal(idx.r, ':')
        is2D = 1;
        % Fixed r: (θ, z) plane
        [Z, Th] = meshgrid(z, th);
        X = Th;
        Y = Z;
        x_label = 'θ [rad]';
        y_label = 'z [m]';
        
    elseif ~isequal(idx.th, ':')
        is2D = 1;
        % Fixed θ: (r, z) plane

        % % Dimensional r:
        % r_frac = linspace(0, 1, par.N_r);  % 0 = inner wall, 1 = outer wall
        % [Z, r_frac_grid] = meshgrid(z, r_frac);  % r_frac_grid is like 0 to 1 vertically
        % 
        % % Now build R from Z
        % R = noz.r_i(Z) + r_frac_grid .* (noz.r_o(Z) - noz.r_i(Z));
        % X = R;

        % Dimensionless r:
        r_frac = linspace(0, 1, par.N_r);  % (0 = inner wall, 1 = outer wall)
        [Z, R_frac] = meshgrid(z, r_frac);  % dimensionless R_frac across Z
        X = R_frac;
        Y = Z;
        x_label = 'r* [-]';
        y_label = 'z [m]';

    else
        is2D = 0;
    end
    
    if is2D == 1
        figure()
        surf(X, Y, T_slice')
        shading interp
        xlabel(x_label)
        ylabel(y_label)
        set(gca, 'YDir', 'reverse')
        set(gca, 'XAxisLocation', 'top')
    
        Tmin = 50;
        Tmax = 1800;
        clim([Tmin, Tmax])
        cb = colorbar;
        ylabel(cb, 'Temperature [K]')
        view([0, 90])

    elseif is2D == 0
        [Zc, Rc, Thc] = ndgrid(z, r, th);
    
        % Densify in z, r, and θ
        z_dense = linspace(min(z),   max(z),   numel(z)*2);   % 4× more z samples
        r_dense = linspace(min(r),   max(r),   numel(r)*1);   % 1× more r samples
        th_dense = linspace(min(th),  max(th),  numel(th)*20);  % 20× more θ samples
        [Zg, Rg, Thg] = ndgrid(z_dense, r_dense, th_dense);
    
        % Interp the temp
        Tg = interpn(Zc, Rc, Thc, T_soln, Zg, Rg, Thg, 'linear');
    
        % Convert to Cartesian (with θ stretch)
        th_scale = 25;
        Xg = Rg .* cos(Thg * th_scale);
        Yg = Rg .* sin(Thg * th_scale);
        Zg = Zg;  % already in Z
    
        % Flatten for scatter3
        Xv = Xg(:);
        Yv = Yg(:);
        Zv = Zg(:);
        Tv = Tg(:);
    
        % Plot
        figure;
        scatter3(Xv, Yv, Zv, 8, Tv, 'filled');
        colormap(jet);  % jet coloring
        colorbar;
        xlabel('x'); ylabel('y'); zlabel('z');
        title('3D Temperature Distribution (Full 3D Interpolation)');
        axis equal tight;
        set(gca, 'ZDir', 'reverse');
        view(3);
    
        % Lighting
        camlight headlight;
        lighting gouraud;
    end
end