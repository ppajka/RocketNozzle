# RocketNozzle

# Instructions for using the MATLAB script:
*** Open the file called MAIN_Rocket_Nozzle_Project.m. This is the main
script — all other functions are auxillary. ***

The script can solve spatially 2D or 3D, and temporally transient or
steady state, so it is important to be accustomed with how to do 
the selections

# Step 1: Number of cells
The user is first asked to provide the number of cells for each
direction in vector form. The input should be given in the following
form: [z, r, θ]

Example:
An input of [10, 5, 3] indicates 10 cells in z, 5 in r, and 3 in θ.

The script runs fairly efficiently, but slows quite a bit for larger
domains. For a run time less than 2 minutes, it is recommended to keep to
the following contraints:

For transient solutions, N_z*N_r*N_th < 15625
For steady-state solutions, N_z*N_r*N_th < 3375
**(system used an i9-13900kf @ ~5.3 GHz, RTX 4070 Ti, and 64 GB of RAM)

In general, transient solution will be much faster (even for large time
spans) because lsqnonlin() is implemented with a naive initial guess and
parameters for the steady-state approach.

Note: if the user wants to perform 2D analysis, then the fixed direction
does not count towards the total cell count, though it still needs to be
specified

# Step 2: Indices of interst
The user is then asked what the indices of interest are. These should
also be given in vector form, [z, r, θ]. The delimiter, ':', is passed in
to indicate that all indices are used for a particular variable. If all
three directions have ':', then 3D analysis is performed on the specified
cells. If the user wants to speed up the analysis, a 2D solution can be
found by fixing one of the variables. For more accurate results, this 2D
slice should be selected such that it intersects the internal channel.
This way, the convective BC propagates information to the interior.

Example:
An input of [5, ':', ':'] indicates a 2D solution at a fixed index value 
of z = 5. An input of [':', ':', ':'] indicates a full 3D solution (with
a longer run time).

# Step 3: Time dependency
The user is prompted if they want a steady-state or transient solution.
Inputting 0 indicates a steady-state using lsqnonlin(), and inputting 1
indicates a transient solution using ode15s().

If transient is selected, the user is prompted for the time span to solve
over.

# Step 4: Solver runs ...

% Step 5: Plotting
By default, the script plots either the steady-state solution, or the
transient solution at the last time step. If a 2D analysis was performed,
a 2D surf plot is automatically generated. If a 3D analysis was 
performed,the user is asked if they want a 3D scatterplot, or if they 
want to specify a 2D slice to plot (using the same convention as step 2).

If other plots are desired, the solution can be saved as a .mat file in
order to be repeatedly used.

