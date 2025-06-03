%% MEGN571: Semester Project 
% Author: Peter Pajka
% Regeneratively cooled rocket nozzle — transient/steady-state 3D/2D

%% Main Script
clc; clear; close all

%% Nozzle Properties
noz.rho = 9130;  % density [kg/m^3]
noz.cp = 385;  % specfic heat, c_p [J/kg-K]
noz.k = 400;  % thermal conductivity [W/m-K]

%% Nozzle Geometric Dimensions
noz.z_start = 1;  % starting length for thermally devloping channel

noz.th_span = pi/50;  % span of nozzle angular sector (100 sectors) [rad]

noz.r_throat = 0.539/2;  % radius of nozzle throat
noz.r_span = 0.04;  % span of nozzle radial thickness [m]
noz.r_i = @(z) noz.r_throat + 1*(z-noz.z_start);  % inner radius (linear profile in z) [m]
noz.r_o = @(z) noz.r_i(z) + noz.r_span;  % outer radius [m]

% Use fsolve to find the nozzle length to achieve 3X throat area
options = optimset('Display', 'off');
noz.z_end = fsolve(@(z) pi*noz.r_i(z)^2 - 3*pi*noz.r_throat^2, 2, options);  % ending length of nozzle [m]
noz.z_span = noz.z_end - noz.z_start;

% Cross-sectional area for combustion gases
noz.A_exh_func = @(z) pi*(noz.r_i(z)).^2;

%% Fluid Properties
flu.rho = 21.97;  % density [kg/m^3]
flu.cp = 15000;  % specfic heat, c_p [J/kg-K]
flu.k = 0.2575;  % thermal conductivity [W/m-K]
flu.mu = 1.16e-5;  % viscosity [Pa-s]

flu.m_dot = 57.166/100;  % mass flow split to 100 channels [kg/s]
flu.T_in = 50;  % entry temp [K]

%% Channel Geometry
chnl.r_span = 0.02;  % span of channel radial thickness [m]
chnl.th_span = pi/100;  % span of channel angular sector [rad]

chnl.r_i = @(z) noz.r_i(z) + 0.005;  % offset length from inner nozzle wall [m]
chnl.r_o = @(z) chnl.r_i(z) + chnl.r_span;  % offset length from outer nozzle wall [m]

% Center the channel in the nozzle:
chnl.th_start = noz.th_span/2 - chnl.th_span/2;  % starting angle for channel
chnl.th_end = noz.th_span/2 + chnl.th_span/2;  % ending angle for channel

% Cross-sectional area of fluid flow
chnl.A_flu_func = @(z) 0.5*chnl.th_span.*(chnl.r_o(z).^2 - chnl.r_i(z).^2);

% Hydraulic Diameter
chnl.D_h = @(z) 2*(chnl.r_o(z).^2 - chnl.r_i(z).^2)*chnl.th_span ...
    ./ (2*(chnl.r_o(z) - chnl.r_i(z)) + chnl.th_span*(chnl.r_o(z) + chnl.r_i(z)));

% Reynolnds Number
chnl.Re = @(z) flu.m_dot*chnl.D_h(z) ./ (flu.mu*chnl.A_flu_func(z));

%% Convection Properties
air.h_T = 25;  % outside convection coeff [W/m^2-K]
air.T = 300;  % outside temp [K]

exh.h_T = 50;  % exhaust convection coeff [W/m^2-K] (placeholder, use Bartz later)
exh.T = 2000;  % inside temp [K] (placeholder)

% Channel convection coeff [W/m^2-K]
flu.h_T_func = @(z) (flu.k ./ z)*0.023 .* (flu.m_dot*chnl.D_h(z) ./ (chnl.A_flu_func(z)*flu.mu)).^0.8 ...
    * (flu.rho*flu.cp/flu.k)^0.4;
flu.T = 500;  % channel temp [K] (placeholder, this is found later)

%% Parameters
Ncell_input = input("Enter the number of cells for z, r, and θ as a vector (e.g., [20, 20, 10]): ");
par.N_z = Ncell_input(1);  % number of cells in z
par.N_r = Ncell_input(2);  % number of cells in r
par.N_th = Ncell_input(3);  % number of cells in θ

% Arrays of z, r, θ values
z = linspace(noz.z_start, noz.z_end, par.N_z);
r = linspace(noz.r_i(z(1)), noz.r_o(z(1)), par.N_r);  % placeholder for mask
th = linspace(0, noz.th_span, par.N_th);

% Boundary axial and radial locations (used for areas)
z_b = [z(1), 0.5*(z(1:end-1) + z(2:end)), z(end)];
r_b = [r(1), 0.5*(r(1:end-1) + r(2:end)), r(end)];  % placeholder
th_b = [th(1), 0.5*(th(1:end-1) + th(2:end)), th(end)];

% Differential spacings for volume/areas
par.dz = z_b(2:end) - z_b(1:end-1);  % spacing in z
par.dr = r_b(2:end) - r_b(1:end-1);  % spacing in r
par.dth = th(2) - th(1);  % spacing in θ


%% Creation of Logical Mask for Channel Geometry
% Find index of cell with boundary closest to channel radii
[~, chnl.r_i_idx] = min(abs(r_b - chnl.r_i(z(1))));
[~, chnl.r_o_idx] = min(abs(r_b - chnl.r_o(z(1))));
chnl.r_o_idx = chnl.r_o_idx - 1;  % outer boundary must shift one index down

% Find index of cell with boundary closest to channel angle
[~, chnl.th_i_idx] = min(abs(th_b - chnl.th_start));
[~, chnl.th_o_idx] = min(abs(th_b - chnl.th_end));
chnl.th_o_idx = chnl.th_o_idx - 1;  % outer boundary must shift one index down

% 3D logical array for masking (1's stay, and 0's are omitted from domain)
mask.mask = true(par.N_z, par.N_r, par.N_th);
mask.mask(:, chnl.r_i_idx:chnl.r_o_idx, chnl.th_i_idx:chnl.th_o_idx) = false;


%% Volume and Areas Calculations
% Pre-allocate radial, axial, and circumferential area matrices
par.A_r = zeros(par.N_z, par.N_r + 1);
par.A_z = zeros(par.N_z + 1, par.N_r);
par.A_th = zeros(par.N_z, par.N_r);

par.rdth = zeros(par.N_z, par.N_r);  % differential changes in θ (r*dθ)
par.V_cell = zeros(par.N_z, par.N_r);  % cell volumes

% Loop over all z locations to find radial/circumferential areas and cell volumes
for i=1:par.N_z
    r = linspace(noz.r_i(0.5*(z_b(i) + z_b(i+1))), noz.r_o(0.5*(z_b(i) + z_b(i+1))), par.N_r);  % array of r values at particular z
        % NOTE: using average between z boundaries to handle first and last half-cells
    r_b = [r(1), 0.5*(r(1:end-1) + r(2:end)), r(end)];  % boundary radii

    par.A_r(i, :) = r_b(1:end)*par.dth*par.dz(i);  % radial areas at boundaries
    par.A_th(i, :) = par.dz(i)*par.dr;  % circumferential areas

    par.V_cell(i, :) = (r_b(2:end).^2 - r_b(1:end-1).^2)*par.dth*par.dz(i);  % cell volumes

    par.rdth(i, :) = 0.5*(r_b(2:end) + r_b(1:end-1))*par.dth;  % differential changes in θ
end

% Loop over all z boundary locations to find axial areas
for i=1:par.N_z+1
    r = linspace(noz.r_i(z_b(i)), noz.r_o(z_b(i)), par.N_r);  % array of r values at particular z boundary
    r_b = [r(1), 0.5*(r(1:end-1) + r(2:end)), r(end)];  % boundary radii

    par.A_z(i, :) = 0.5*(r_b(2:end).^2 - r_b(1:end-1).^2)*par.dth;  % axial areas at boundaries
end


%% Solve
% Define mask indices for initial condtions and solution output (only needed for 2D testing)
idx_input = input("Enter the indices of interest for z, r, and θ as a vector (e.g., 2D: [':', ':', 6], or 3D: [':', ':', ':']): ");

% Store the values in the idx struct
idx.z = idx_input(1);
idx.r = idx_input(2);
idx.th = idx_input(3);

% Mask calcs used for initial guess
chunk = mask.mask(idx.z, idx.r, idx.th);
mask.vec = chunk(:);
mask.idx = find(mask.vec);

% Initial condition (3D array)
T_0 = 350*ones(par.N_z, par.N_r, par.N_th);  % initial guess for strip of θ
T_0_mask = T_0.*mask.mask(idx.z, idx.r, idx.th);  % overlay mask

% Initial guess for fluid
T_0_flu = flu.T_in*ones(par.N_z, 1);

% Extract initial guess vector
chunk = T_0_mask(idx.z, idx.r, idx.th);
T_0_vec = chunk(:);  % convert 3D array to vector
T_0_vec = T_0_vec(mask.idx);  % convert full vector to compact form
T_0_vec = [T_0_vec; T_0_flu];

% Pre-allocate 3D array for soln
T_soln = zeros(par.N_z, par.N_r, par.N_th);

soln_type = input("Transient or steady-state solution (0 = steady, 1 = transient): ");
if soln_type == 0
    disp("Solving...")

    % Steady-state function
    noz_ss = @(T) nozzle_dTdt(z, 0, T, noz, flu, chnl, air, exh, par, mask, idx);
    
    % Solve using lsqnonlin
    opts = optimoptions('lsqnonlin', 'Display', 'iter', 'MaxFunctionEvaluations', 1e5);
    T_soln_vec = lsqnonlin(noz_ss, T_0_vec, [], [], opts);
    
    flu.T_soln = T_soln_vec(end - par.N_z + 1:end);  % fluid T is at end of vector
    T_soln_vec = T_soln_vec(1:end - par.N_z);  % nozzle T is at beginning of vector

elseif soln_type == 1
    % Time span for soln
    t_max = input("Enter the time duration (s): ");
    t_span = [0, t_max];  % time span
    disp("Solving...")
    
    % Transient function
    noz_trans = @(t, T) nozzle_dTdt(z, t, T, noz, flu, chnl, air, exh, par, mask, idx);

    % ODE solver
    [t, T_soln_vec] = ode15s(noz_trans, t_span, T_0_vec);
    
    % Extract channel and nozzle temps from T vector
    T_soln_vec = T_soln_vec(end, :);  % grab temp field at final time
    flu.T_soln = T_soln_vec(end - par.N_z + 1:end);  % fluid T is at end of vector
    T_soln_vec = T_soln_vec(1:end - par.N_z);  % nozzle T is at beginning of vector

end

% Insert T vector into full 3D array
T_soln_vec_full(mask.idx) = T_soln_vec;  % convert compact vector to full form
chunk = T_soln(idx.z, idx.r, idx.th);  % chunk extracted for particular slice

% Solution as 3D array:
T_soln(idx.z, idx.r, idx.th) = reshape(T_soln_vec_full(:), size(chunk));

% Insert T_m of fluid into full solution
flu.T_fill = repmat(flu.T_soln, 1, par.N_r, par.N_th);
T_soln(mask.mask == 0) = flu.T_fill(mask.mask == 0);

%% Plotting
n_fixed = (~isequal(idx.z, ':')) + (~isequal(idx.r, ':')) + (~isequal(idx.th, ':'));
is2D = (n_fixed == 1);  % check if any index input is a scalar (indicating 2D)

if is2D == 0  % 3D analysis
    idx_input = input("Enter the indices for plotting (e.g., 3D scatter plot: [':', ':', ':'],  2D surf plot: [8, ':', ':']): ");
    % Store the values in the idx struct
    idx.z = idx_input(1);
    idx.r = idx_input(2);
    idx.th = idx_input(3);
end

plot_temp(T_soln, noz, par, z, r, th, idx)

%% Function for calculating dT/dt
function [dT_dt_vec] = nozzle_dTdt(z, t, T_vec, noz, flu, chnl, air, exh, par, mask, idx)
    % --- Extract parameters ---
    % Number of cells
    N_z = par.N_z;
    N_r = par.N_r;
    N_th = par.N_th;
    % Cell spacings
    dz = par.dz(2);  % use full cell diff (not half-cell)
    dr = par.dr(2);
    rdth = par.rdth;  % matrix rdth(z, r)
    % Areas and Volumes
    A_z = par.A_z;
    A_r = par.A_r;
    A_th = par.A_th;
    V_cell = par.V_cell;
    

    % --- Extract channel and nozzle temps from T vector ---
    flu.T = T_vec(end - N_z + 1:end);  % fluid T is at end of vector
    noz.T_vec = T_vec(1:end - N_z);  % nozzle T is at beginning of vector

    % --- Insert nozzle T vector into full 3D array ---
    noz.T_vec_full = zeros(size(mask.vec));  % pre-allocate full temp vector
    noz.T_vec_full(mask.idx) = noz.T_vec;  % convert compact vector into full form
    T = zeros(N_z, N_r, N_th);  % pre-allocate T (3D array)
    chunk = T(idx.z, idx.r, idx.th);  % chunk being extracted for discretization
    T(idx.z, idx.r, idx.th) = reshape(noz.T_vec_full, size(chunk));  % convert vector to 3D array

    % --- Logic for deciding 3D or 2D analysis ---
    % ***This section of code duplicates a 2D slice such that conduction is
    % zero in the third direction. This only occurs if one index is
    % fixed.
    n_fixed = (~isequal(idx.z, ':')) + (~isequal(idx.r, ':')) + (~isequal(idx.th, ':'));
    is2D = (n_fixed == 1);  % check if any index input is a scalar (indicating 2D)

    if is2D
        if ~isequal(idx.z, ':')
            % extract the 2D slice that was computed
            T_slice = T(idx.z, :, :);
            % replicate along third dimension
            T = repmat(T_slice, [N_z, 1, 1]);
        elseif ~isequal(idx.r, ':') 
            T_slice = T(:, idx.r, :);
            T = repmat(T_slice, [1, N_r, 1]);
        elseif ~isequal(idx.th, ':')
            T_slice = T(:, :, idx.th);
            T = repmat(T_slice, [1, 1, N_th]);
        end
    end

    % --- Find exhaust heat transfer coeff and temp based on initial condition ---
    noz.T_w = mean(T(:, 1, :), 3)';  % LHS wall temperature averaged over all θ (row vec)
    noz.A_exh = noz.A_exh_func(z);  % X-section areas (row vec)
    [exh.h_T, exh.T] = Bartz(noz.T_w, noz.A_exh);
    exh.T = exh.T';  % exhaust temp (col vec)
    exh.h_T = exh.h_T';  % exhaust heat trans coeff (col vec for discretization)

    % --- Find channel mean wall temp and fluid heat transfer coefficient ---
    % Edges on border of channel
    T_edges = {
        T(:, chnl.r_i_idx-1, chnl.th_i_idx:chnl.th_o_idx);  % edge at r=r_i
        T(:, chnl.r_o_idx+1, chnl.th_i_idx:chnl.th_o_idx);  % edge at r=r_o
        T(:, chnl.r_i_idx:chnl.r_o_idx, chnl.th_i_idx-1);  % edge at th=th_i
        T(:, chnl.r_i_idx:chnl.r_o_idx, chnl.th_o_idx+1);  % edge at th=th_o
    };

    T_valid_edges = {};  % valid edges (in the case of 2D analysis there are erroneous zeros)
    for i = 1:4
        T_edge = reshape(T_edges{i}, size(T, 1), []);  % flatten into (Nz, num_points)
        if any(abs(T_edge(:)) > 1e-12)  % check if this edge has any real values
            T_valid_edges{end+1} = T_edge;  % keep it
        end
    end
    
    T_valid_edges = cat(2, T_valid_edges{:});  % concatenate all valid edges
    noz.T_m = mean(T_valid_edges, 2);  % average to get T_m of nozzle
    flu.h_T = flu.h_T_func(z)';  % fluid heat trans coeff (col vec)


    % --- 3D dicretization ---
    % Initialize dT/dt array
    dT_dt = zeros(N_z, N_r, N_th);
    
    % --- Edge boundary conditions ---
    % BC: z=1, r=2:N_r-1 (adiabatic)
    BC_idx = {1, 2:N_r-1, 1:N_th};
    dT_dt(BC_idx{1}, BC_idx{2}, BC_idx{3}) = ...
        (1 ./ (noz.rho*noz.cp*V_cell(BC_idx{1}, BC_idx{2}))) .* ...
        (0 ... %%% BC
        + cond(par, noz.k, T, A_z, BC_idx, dz, 1, 1) ...
        + cond(par, noz.k, T, A_r, BC_idx, dr, 2, 0) ...
        + cond(par, noz.k, T, A_r, BC_idx, dr, 2, 1) ...
        + cond(par, noz.k, T, A_th, BC_idx, rdth, 3, 0) ...
        + cond(par, noz.k, T, A_th, BC_idx, rdth, 3, 1));

    % BC: z=N_z, r=2:N_r-1 (adiabatic)
    BC_idx = {N_z, 2:N_r-1, 1:N_th};
    dT_dt(BC_idx{1}, BC_idx{2}, BC_idx{3}) = ...
        (1 ./ (noz.rho*noz.cp*V_cell(BC_idx{1}, BC_idx{2}))) .* ...
        (cond(par, noz.k, T, A_z, BC_idx, dz, 1, 0) ...
        + 0 ... %%% BC
        + cond(par, noz.k, T, A_r, BC_idx, dr, 2, 0) ...
        + cond(par, noz.k, T, A_r, BC_idx, dr, 2, 1) ...
        + cond(par, noz.k, T, A_th, BC_idx, rdth, 3, 0) ...
        + cond(par, noz.k, T, A_th, BC_idx, rdth, 3, 1));

    % BC: z=2:N_z-1, r=1 (exhaust convection)
    BC_idx = {2:N_z-1, 1, 1:N_th};
    dT_dt(BC_idx{1}, BC_idx{2}, BC_idx{3}) = ...
        (1 ./ (noz.rho*noz.cp*V_cell(BC_idx{1}, BC_idx{2}))) .* ...
        (cond(par, noz.k, T, A_z, BC_idx, dz, 1, 0) ...
        + cond(par, noz.k, T, A_z, BC_idx, dz, 1, 1) ...
        + conv(exh.h_T(BC_idx{1}), T, A_r, BC_idx, exh.T(BC_idx{1}), 2, 0) ... %%% BC
        + cond(par, noz.k, T, A_r, BC_idx, dr, 2, 1) ...
        + cond(par, noz.k, T, A_th, BC_idx, rdth, 3, 0) ...
        + cond(par, noz.k, T, A_th, BC_idx, rdth, 3, 1));

    % BC: z=2:N_z-1, r=N_r (air convection)
    BC_idx = {2:N_z-1, N_r, 1:N_th};
    dT_dt(BC_idx{1}, BC_idx{2}, BC_idx{3}) = ...
        (1 ./ (noz.rho*noz.cp*V_cell(BC_idx{1}, BC_idx{2}))) .* ...
        (cond(par, noz.k, T, A_z, BC_idx, dz, 1, 0) ...
        + cond(par, noz.k, T, A_z, BC_idx, dz, 1, 1) ...
        + cond(par, noz.k, T, A_r, BC_idx, dr, 2, 0) ...
        + conv(air.h_T, T, A_r, BC_idx, air.T, 2, 1) ... %%% BC
        + cond(par, noz.k, T, A_th, BC_idx, rdth, 3, 0) ...
        + cond(par, noz.k, T, A_th, BC_idx, rdth, 3, 1));


    % --- Corner boundary conditions ---
    % BC: z=1, r=1 (adiabatic + exhaust convection)
    BC_idx = {1, 1, 1:N_th};
    dT_dt(BC_idx{1}, BC_idx{2}, BC_idx{3}) = ...
        (1 ./ (noz.rho*noz.cp*V_cell(BC_idx{1}, BC_idx{2}))) .* ...
        (0 ... %%% BC
        + cond(par, noz.k, T, A_z, BC_idx, dz, 1, 1) ...
        + conv(exh.h_T(BC_idx{1}), T, A_r, BC_idx, exh.T(BC_idx{1}), 2, 0) ... %%% BC
        + cond(par, noz.k, T, A_r, BC_idx, dr, 2, 1) ...
        + cond(par, noz.k, T, A_th, BC_idx, rdth, 3, 0) ...
        + cond(par, noz.k, T, A_th, BC_idx, rdth, 3, 1));

    % BC: z=1, r=N_r (adiabatic + air convection)
    BC_idx = {1, N_r, 1:N_th};
    dT_dt(BC_idx{1}, BC_idx{2}, BC_idx{3}) = ...
        (1 ./ (noz.rho*noz.cp*V_cell(BC_idx{1}, BC_idx{2}))) .* ...
        (0 ... %%% BC
        + cond(par, noz.k, T, A_z, BC_idx, dz, 1, 1) ...
        + cond(par, noz.k, T, A_r, BC_idx, dr, 2, 0) ...
        + conv(air.h_T, T, A_r, BC_idx, air.T, 2, 1) ... %%% BC
        + cond(par, noz.k, T, A_th, BC_idx, rdth, 3, 0) ...
        + cond(par, noz.k, T, A_th, BC_idx, rdth, 3, 1));

    % BC: z=N_z, r=1 (adiabatic + exhaust convection)
    BC_idx = {N_z, 1, 1:N_th};
    dT_dt(BC_idx{1}, BC_idx{2}, BC_idx{3}) = ...
        (1 ./ (noz.rho*noz.cp*V_cell(BC_idx{1}, BC_idx{2}))) .* ...
        (cond(par, noz.k, T, A_z, BC_idx, dz, 1, 0) ...
        + 0 ... %%% BC
        + conv(exh.h_T(BC_idx{1}), T, A_r, BC_idx, exh.T(BC_idx{1}), 2, 0) ... %%% BC
        + cond(par, noz.k, T, A_r, BC_idx, dr, 2, 1) ...
        + cond(par, noz.k, T, A_th, BC_idx, rdth, 3, 0) ...
        + cond(par, noz.k, T, A_th, BC_idx, rdth, 3, 1));

    % BC: z=N_z, r=N_r (adiabatic + air convection)
    BC_idx = {N_z, N_r, 1:N_th};
    dT_dt(BC_idx{1}, BC_idx{2}, BC_idx{3}) = ...
        (1 ./ (noz.rho*noz.cp*V_cell(BC_idx{1}, BC_idx{2}))) .* ...
        (cond(par, noz.k, T, A_z, BC_idx, dz, 1, 0) ...
        + 0 ... %%% BC
        + cond(par, noz.k, T, A_r, BC_idx, dr, 2, 0) ...
        + conv(air.h_T, T, A_r, BC_idx, air.T, 2, 1) ... %%% BC
        + cond(par, noz.k, T, A_th, BC_idx, rdth, 3, 0) ...
        + cond(par, noz.k, T, A_th, BC_idx, rdth, 3, 1));


    % --- Interior cells ---
    int_idx = {2:N_z-1, 2:N_r-1, 1:N_th};
    dT_dt(int_idx{1}, int_idx{2}, int_idx{3}) = ...
        (1 ./ (noz.rho*noz.cp*V_cell(int_idx{1}, int_idx{2}))) .* ...
        (cond(par, noz.k, T, A_z, int_idx, dz, 1, 0) ...
        + cond(par, noz.k, T, A_z, int_idx, dz, 1, 1) ...
        + cond(par, noz.k, T, A_r, int_idx, dr, 2, 0) ...
        + cond(par, noz.k, T, A_r, int_idx, dr, 2, 1) ...
        + cond(par, noz.k, T, A_th, int_idx, rdth, 3, 0) ...
        + cond(par, noz.k, T, A_th, int_idx, rdth, 3, 1));


    % --- Overwrite cells bordering masked cells ---
    % BC: r=chnl.r_i (convection to channel)
    chnl_idx = {1:N_z, chnl.r_i_idx-1, chnl.th_i_idx:chnl.th_o_idx};  % cells before omission
    dT_dt(chnl_idx{1}, chnl_idx{2}, chnl_idx{3}) = ...
        dT_dt(chnl_idx{1}, chnl_idx{2}, chnl_idx{3}) ...
        + (1 ./ (noz.rho*noz.cp*V_cell(chnl_idx{1}, chnl_idx{2}))) .* ...
        (-cond(par, noz.k, T, A_r, chnl_idx, dr, 2, 1) ...  % subtract off already set conduction
        + conv(flu.h_T(chnl_idx{1}), T, A_r, chnl_idx, flu.T(chnl_idx{1}), 2, 1));  % add on convection to replace conduction

    % BC: r=chnl.r_o (convection to channel)
    chnl_idx = {1:N_z, chnl.r_o_idx+1, chnl.th_i_idx:chnl.th_o_idx};  % cells before omission
    dT_dt(chnl_idx{1}, chnl_idx{2}, chnl_idx{3}) = ...
        dT_dt(chnl_idx{1}, chnl_idx{2}, chnl_idx{3}) ...
        + (1 ./ (noz.rho*noz.cp*V_cell(chnl_idx{1}, chnl_idx{2}))) .* ...
        (-cond(par, noz.k, T, A_r, chnl_idx, dr, 2, 0) ...
        + conv(flu.h_T(chnl_idx{1}), T, A_r, chnl_idx, flu.T(chnl_idx{1}), 2, 0));

    % BC: th=chnl.th_i (convection to channel)
    chnl_idx = {1:N_z, chnl.r_i_idx:chnl.r_o_idx, chnl.th_i_idx-1};  % cells before omission
    dT_dt(chnl_idx{1}, chnl_idx{2}, chnl_idx{3}) = ...
        dT_dt(chnl_idx{1}, chnl_idx{2}, chnl_idx{3}) ...
        + (1 ./ (noz.rho*noz.cp*V_cell(chnl_idx{1}, chnl_idx{2}))) .* ...
        (-cond(par, noz.k, T, A_th, chnl_idx, rdth, 3, 1) ...
        + conv(flu.h_T(chnl_idx{1}), T, A_th, chnl_idx, flu.T(chnl_idx{1}), 3, 1));

    % BC: th=chnl.th_o (convection to channel)
    chnl_idx = {1:N_z, chnl.r_i_idx:chnl.r_o_idx, chnl.th_o_idx+1};  % cells before omission
    dT_dt(chnl_idx{1}, chnl_idx{2}, chnl_idx{3}) = ...
        dT_dt(chnl_idx{1}, chnl_idx{2}, chnl_idx{3}) ...
        + (1 ./ (noz.rho*noz.cp*V_cell(chnl_idx{1}, chnl_idx{2}))) .* ...
        (-cond(par, noz.k, T, A_th, chnl_idx, rdth, 3, 0) ...
        + conv(flu.h_T(chnl_idx{1}), T, A_th, chnl_idx, flu.T(chnl_idx{1}), 3, 0));


    % --- Fluid mean temp calcs ---
    flu.dT_dt = zeros(N_z, 1);

    % --- Fluid initial condition ---
    % BC: z=1
    flu_idx = 1;
    flu.dT_dt(flu_idx) = flu_mean_temp(z, noz, flu_idx, flu, chnl, par);

    % --- Fluid interior cells ---
    flu_idx = 2:N_z;
    flu.dT_dt(flu_idx) = flu_mean_temp(z, noz, flu_idx, flu, chnl, par);


    % --- Extract T vector from 3D array ---
    dT_dt_mask = dT_dt.*mask.mask(idx.z, idx.r, idx.th);  % mask over dT_dt
    chunk = dT_dt_mask(idx.z, idx.r, idx.th);  % chunk being extracted for ODE solver
    dT_dt_vec = chunk(:);  % convert 3D array to vector
    dT_dt_vec = dT_dt_vec(mask.idx);  % convert full vector to compact form
    dT_dt_vec = [dT_dt_vec; flu.dT_dt];  % append fluid dT/dt vector to end
end
