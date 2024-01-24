function Z_I_star = initialise_approximation_structure(params, fem, I_star)
%INITIALISE_APPROXIMATION_STRUCTURE Initialises a structure containing
%collocation points, and associated FEM and timestepping information.
% Input:
%   params     Structure with approximation parameters
%   fem        Structure with finite element data
%   I_star     Structure with collocation points in enhanced grid.
% Output: 
%   Z_I_star   Structure with grid, approximation and times

    % Base structure on reduced grid
    I_star_r = reduce_sparse_grid(I_star);
    Z_I_star = cell(I_star_r.size,1);
    % For each collocation point initialise data.
    for ii = 1:I_star_r.size
        z = I_star_r.knots(:,ii);
        z_struct.z = z;
        
        % Timestepping data
        z_struct.u_z = zeros(size(fem.xy,1),1);
        z_struct.u_z_tplusdt = z_struct.u_z;
        z_struct.dt_z_tplusdt = params.dt0;
        z_struct.t_z = 0;
        z_struct.tplusdt = 0;

        % Advection-diffusion fields
        z_struct.a_fn = @(x1,x2) fem.a_fn(x1,x2,z);
        z_struct.wind_fn = @(x1,x2) fem.wind_fn(x1,x2,z);

        % FEM matrices and vectors
        z_struct.Q = fem.mass;
        z_struct.D = fem.diff(z);
        z_struct.N = fem.conv(z);
        z_struct.f = fem.f; %@(t) fem.f(t,z);
        z_struct.bc = @(t) fem.bc_fn(t,z);
        z_struct.bc_prime = @(t) fem.bc_fn_prime(t,z);
        z_struct.bound = fem.bound;
        z_struct.notbound = fem.notbound;

        z_struct.dt_z = params.dt0;
        z_struct.w = nan(length(fem.bound),1);
        z_struct.n_steps = 0;
        z_struct.ge_estimate = 0;

        Z_I_star{ii} = z_struct;
    end
    
end
