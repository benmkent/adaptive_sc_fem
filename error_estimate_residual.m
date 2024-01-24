function [error, precompute] = error_estimate_residual(error, t_r, delta_t, Z_I_star, I_star, I, params, fem, precompute)
%ERROR_ESTIMATE_RESIDUAL Computes the residual error estimators
%
% Inputs        error           Previous error structure
%               t_r             Current time
%               delta_t         Algorithm timestep
%               Z_I_star        data structure of collocation points
%               I_star          Enhanced sparse grid
%               I               Sparse grid
%               params          Approximation parameters
%               fem             Finite element data structure
%               precompute      Precomputed data
% Outputs       error           Structure of computed error estimates
%               Z_I_star        Updated data structure

%% Initialise variables
% Find I in I star
I_r = reduce_sparse_grid(I);
I_star_r = reduce_sparse_grid(I_star);
if I_star_r.size == I_r.size
    I_in_I_star = 1:I_r.size;
else
    [~,~,I_in_I_star] = compare_sparse_grids(I_star, I_star_r, I, I_r);
end

eta_x_z=[];
eta_x_z_T=[];
eta_u_z = [];
eta_t_z=[];
eta_t_k_z=[];
eta_t_z_du=[];
eta_t_z_u_w=[];
eta_t_z_eta_w=[];

%% First compute spatial and temporal terms for each collocation point
% Note that it may be worthwhile to combine and vectorise the spatial and
% temporal error estimators as there is an inner loop over each timestep.
for ii_z = 1:length(I_in_I_star)
    fprintf('...pt %d of %d', ii_z, length(Z_I_star));
    ii_in_I_star = I_in_I_star(ii_z);

    %% First extract only the time intervals between t_r=t and t_{r+1}=t+timestepAlg
    % If set up as in Chapter 5  this should be all intervals.
    % If timesteps are not forced to align to sync times this is a little
    % more complicated
    z_for_error_est = Z_I_star{ii_in_I_star};
    % Find indices between t and t+delta_t
    inds_for_est = (t_r<z_for_error_est.t_z) & (z_for_error_est.t_z < t_r + delta_t);
    % Now get the index before and after
    ind_before_or_equal_to_t = find((t_r<z_for_error_est.t_z),1,'first')-1;
    ind_greater_or_equal_to_tplusdt = find((z_for_error_est.t_z < t_r + delta_t),1,'last')+1;
    % Construct approximation on times t_r,0 , t_r,1,....,t_r,n
    t_prev = z_for_error_est.t_z(ind_before_or_equal_to_t);
    u_prev = z_for_error_est.u_z(:,ind_before_or_equal_to_t);
    t_r1 = z_for_error_est.t_z(ind_before_or_equal_to_t+1);
    u_tr1 = z_for_error_est.u_z(:,ind_before_or_equal_to_t+1);

    t_rnm1 = z_for_error_est.t_z(ind_greater_or_equal_to_tplusdt-1);
    u_rnm1 = z_for_error_est.u_z(:,ind_greater_or_equal_to_tplusdt-1);
    t_next = z_for_error_est.t_z(ind_greater_or_equal_to_tplusdt);
    u_next = z_for_error_est.u_z(:,ind_greater_or_equal_to_tplusdt);

    % If not aligned with sync times we must pass the full steps and then
    % integrate the first and final steps over the appropriate
    % fraction of the interval.
    z_for_error_est.u_z = ...
        [u_prev, z_for_error_est.u_z(:,inds_for_est), u_next];
    z_for_error_est.t_z = [t_prev; z_for_error_est.t_z(inds_for_est); t_next];
    z_for_error_est.dt_z = [z_for_error_est.dt_z(ind_before_or_equal_to_t);...
        z_for_error_est.dt_z(inds_for_est);...
        z_for_error_est.dt_z(ind_greater_or_equal_to_tplusdt)];

    %% Spatial error estimators
    % Compute spatial estimates
    [eta_x_z(ii_z),eta_x_z_T{ii_z}]  = compute_element_spatial_errors(z_for_error_est, fem,t_r,t_r+delta_t);
    [eta_u_z{ii_z}] = Z_I_star{ii_in_I_star}.u_z_tplusdt;
    fprintf('...x done')

    %% Temporal error estimators
    % Compute temporal terms
    [eta_t_z(ii_z), eta_t_k_z{ii_z}, eta_t_z_du{ii_z},...
        eta_t_z_u_w{ii_z}, eta_t_z_eta_w{ii_z}] = compute_temporal_errors(z_for_error_est, fem,t_r,t_r+delta_t);
    fprintf('...t done\n');
end

%% Interpolation error estimator
% Compute interpolation error estimator
fprintf('Compute interpolation error\n');
[eta_y,eta_y_mi, eta_mi] = compute_interpolation_errors(t_r,t_r+delta_t,Z_I_star, I, I_star, fem, params, precompute);
fprintf('...done\n');

%% Compute complete error estimators
% Define constants
Lebesgue_r = (I_r.size)^2; % Chifka Lemma 3.1
amax = fem.amax;
amin = fem.amin;

pi_x_r = sqrt(1/amin * 6*Lebesgue_r * sum(eta_x_z.^2 .* precompute.l1_norms_I));
pi_t_r = sqrt(1/amin * 6*Lebesgue_r *(1+amax^2) * sum(eta_t_z.^2 .* precompute.l1_norms_I));
pi_y_r = sqrt(1/amin * 2*eta_y^2);

%% Write to output structure
error.pi_x_r = pi_x_r;
error.pi_t_r = pi_t_r;
error.pi_y_r = pi_y_r;
error.pi_y_mi_r = eta_y_mi;
error.y_mi_r = eta_mi;
error.eta_x_z_T = eta_x_z_T;
error.eta_u_z = eta_u_z;
error.eta_t_z_du = eta_t_z_du;
error.eta_t_z_u_w = eta_t_z_u_w;
error.eta_t_z_eta_w = eta_t_z_eta_w;
end
