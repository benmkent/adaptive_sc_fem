%% Generate path
addpath('./subroutines')
addpath('./tools')
addpath('./ifiss_modified')
addpath('./verfurth/')
addpath('./vidlickova/')
addpath('./sparse-grid')
addpath('./bk')

addpath(genpath('../ifiss3.6'))
addpath(genpath('../tifiss1.2'))
addpath(genpath('../sparse-grids-matlab-kit_v-22-02_California'))

%% Define inputs
if ~exist('problem','var')
    problem = define_problem('doubleglazing');
end
if ~exist('params','var')
    params = define_params();
end
if ~exist('reference','var')
    reference = define_reference();
end

%% Initialise
fprintf('Initialise sparse grid...'); tic();
[I_star, I] = initialise_sparse_grid(problem, params);
I_star_r = reduce_sparse_grid(I_star);
fprintf('done (%f seconds)\n', toc())

% FEM approximation
fprintf('Initialise FEM matrices...'); tic();
% fem = initialise_fem_matrices(problem,params);
fem = initialise_fem_matrices_new(problem,params);
fprintf('done (%f seconds)\n', toc());
fprintf('Total DOF %d\n', size(fem.xy,1) * I_star_r.size);

fprintf('Initialise approximation structure...'); tic();
Z_I_star = initialise_approximation_structure(problem, params, fem, I_star);
params_lofi = params; params_lofi.letol = params.letol_lofi; params_lofi.dt0 = params.dt0_lofi;
Z_I_star_lofi = initialise_approximation_structure(problem, params_lofi, fem, I_star);
fprintf('done (%f seconds)\n', toc())

fprintf('Initialise precomputed data (Lagrange product integrals)...'); tic();
precompute = initialise_precomputed_data(I_star, I, params, problem);
fprintf('done (%f seconds)\n', toc())

error_vidlickova = initialise_vidlickova_error();

delta_t = params.t0;
k_delta_t = 2;

%% Initilaise datatable
data_table = write_to_data_table();
%% Initialise figures
axArray = initialise_figures();

%% Loop
t=0; E = 1e-5;
reference_times = reference.times;

piV =0; piV_I=0; piV_delta = 0; piV_space =0; pi_x = 0; pi_t = 0; pi_y = 0;
while t < params.T
    %% SOLVE
    fprintf('Propagate SC point solutions to at least t=%f...\n',t+delta_t); tic();
    for ii_z = 1:length(Z_I_star)
        fprintf('.... (target t= %f) %d of %d',t+delta_t, ii_z, length(Z_I_star));
        Z_I_star{ii_z} = propagate_sc_point( Z_I_star{ii_z}, t + delta_t, params, fem);
        Z_I_star_lofi{ii_z} = propagate_sc_point(Z_I_star_lofi{ii_z}, t + delta_t, params_lofi, fem);
        dt_tplusdt = Z_I_star{ii_z}.dt_z_tplusdt;
        fprintf('...(%f timestep at target t)', dt_tplusdt);
        fprintf('...done\n');
    end
    fprintf('done (%f seconds)\n', toc())

    %% ESTIMATE
    fprintf('Compute error estimates...'); tic();
    fprintf('.... Precompute vector-matrix-vector products\n')
    precompute = precompute_interp_and_l2spatial(Z_I_star, I_star, I, fem, precompute);

    if params.jomp_estimator == 1
        fprintf('Compute BK hierarchical estimators...\n')
        [error_bk, Z_I_star, precompute] = error_estimate_bk(Z_I_star, Z_I_star_lofi, I_star, I, params, fem, precompute);
    else
        error_bk.pi = nan;
        error_bk.pi_I = 0;
        error_bk.pi_delta = inf;
        error_bk.pi_I_delta= inf;
        error_bk.pi_I_alpha = {};
    end
    error{2} = error_bk;

    if params.residual_estimator == 1
        fprintf('Compute spatial and temporal residual estimators\n')
        % Find I in I star
        I_r = reduce_sparse_grid(I);
        I_star_r = reduce_sparse_grid(I_star);
        [~,~,I_in_I_star] = compare_sparse_grids(I_star, I_star_r, I, I_r);
        for ii_z = 1:length(I_in_I_star)
            fprintf('...pt %d of %d', ii_z, length(Z_I_star));
            ii_in_I_star = I_in_I_star(ii_z);

            %%%%%%%%%%%
            % Must first extract only the time intervals between t_r=t and t_{r+1}=t+timestepAlg
            %%%%%%%%%%%
            z_for_error_est = Z_I_star{ii_in_I_star};
            % Find indices between t and t+delta_t
            inds_for_est = (t<z_for_error_est.t_z) & (z_for_error_est.t_z < t + delta_t);
            % Now get the index before and after
            ind_before_or_equal_to_t = find((t<z_for_error_est.t_z),1,'first')-1;
            ind_greater_or_equal_to_tplusdt = find((z_for_error_est.t_z < t + delta_t),1,'last')+1;
            % Construct approximation on times t_r,0 , t_r,1,....,t_r,n
            t_prev = z_for_error_est.t_z(ind_before_or_equal_to_t);
            u_prev = z_for_error_est.u_z(:,ind_before_or_equal_to_t);
            t_r1 = z_for_error_est.t_z(ind_before_or_equal_to_t+1);
            u_tr1 = z_for_error_est.u_z(:,ind_before_or_equal_to_t+1);

            t_rnm1 = z_for_error_est.t_z(ind_greater_or_equal_to_tplusdt-1);
            u_rnm1 = z_for_error_est.u_z(:,ind_greater_or_equal_to_tplusdt-1);
            t_next = z_for_error_est.t_z(ind_greater_or_equal_to_tplusdt);
            u_next = z_for_error_est.u_z(:,ind_greater_or_equal_to_tplusdt);
          
%             z_for_error_est.u_z = ...
%                 [interp1([t_prev,t_r1],[u_prev, u_tr1].', t).',...
%                     z_for_error_est.u_z(:,inds_for_est),...
%                     interp1([t_rnm1, t_next],[u_rnm1, u_next].', t+delta_t).'...
%                     ];
%             z_for_error_est.t_z = [t; z_for_error_est.t_z(inds_for_est); t+delta_t];
            % Don't use interpolation. We must pass the full steps and then
            % integrate the first and final steps over the appropriate
            % fraction of the interval.
            z_for_error_est.u_z = ...
                [u_prev, z_for_error_est.u_z(:,inds_for_est), u_next];
            z_for_error_est.t_z = [t_prev; z_for_error_est.t_z(inds_for_est); t_next];
            z_for_error_est.dt_z = [z_for_error_est.dt_z(ind_before_or_equal_to_t);...
                z_for_error_est.dt_z(inds_for_est);...
                z_for_error_est.dt_z(ind_greater_or_equal_to_tplusdt)];

            %%%%%%%%%%%%
            % Spatial error
            %%%%%%%%%%%%
            % Compute spatial estimates eqn (6.16)
%             [eta_x_z(ii_z), eta_x_k_z{ii_z}, eta_x_z_T{ii_z}, eta_x_k_z_T{ii_z},eta_z_lp(ii_z),eta_z_T_lp{ii_z}] = compute_element_spatial_errors(Z_I_star{ii_in_I_star}, fem);
            [eta_x_z(ii_z),eta_x_z_T{ii_z}]  = compute_element_spatial_errors_new(z_for_error_est, fem,t,t+delta_t);
            [eta_u_z{ii_z}] = Z_I_star{ii_in_I_star}.u_z_tplusdt;
            fprintf('...x done')
            
            %%%%%%%%%%%%%%
            % Temporal error
            %%%%%%%%%%%%%%
            % Compute temporal terms Lemma 4.2 and (4.25)
%             [eta_t_z(ii_z), eta_t_k_z{ii_z}, dU_k_z{ii_z}, u_w_k_z{ii_z}, eta_w_k_z{ii_z}] = compute_temporal_estimator(Z_I_star{ii_in_I_star}, fem);
            [eta_t_z(ii_z), eta_t_k_z{ii_z}, eta_t_z_du{ii_z}, eta_t_z_u_w{ii_z}, eta_t_z_eta_w{ii_z}] = compute_temporal_errors_new(z_for_error_est, fem,t,t+delta_t);
            fprintf('...t done\n');
        end

        %%%%%%%%%%%%%%
        %  Interpolation error
        %%%%%%%%%%%%%%
        fprintf('Compute interpolation error\n');
%         [eta_y_mi] = compute_interpolation_estimator(Z_I_star, I, I_star, fem, params, precompute);
        [eta_y,eta_y_mi, eta_mi] = compute_interpolation_errors(t,t+delta_t,Z_I_star, I, I_star, fem, params, precompute);
        fprintf('...done\n');

        %%%%%%%%%%%%%%
        % Compute error estimators
        %%%%%%%%%%%%%%
        % Compute spatial error eqn (4.11)
        Lebesgue_r = (I_r.size)^2; % Chifka Lemma 3.1  
        amax = 0.1;
        amin = 0.1;

        pi_x_r = sqrt(1/amin * 6*Lebesgue_r * sum(eta_x_z.^2 .* precompute.l1_norms_I));
        pi_t_r = sqrt(1/amin * 6*Lebesgue_r *(1+amax^2) * sum(eta_t_z.^2 .* precompute.l1_norms_I));
        pi_y_r = sqrt(1/amin * 2*eta_y^2);

        error{3}.pi_x_r = pi_x_r;
        error{3}.pi_t_r = pi_t_r;
        error{3}.pi_y_r = pi_y_r;
        error{3}.pi_y_mi_r = eta_y_mi;
        error{3}.y_mi_r = eta_mi;
        error{3}.eta_x_z_T = eta_x_z_T;
        error{3}.eta_u_z = eta_u_z;
        error{3}.eta_t_z_du = eta_t_z_du;
        error{3}.eta_t_z_u_w = eta_t_z_u_w;
        error{3}.eta_t_z_eta_w = eta_t_z_eta_w;

        error{3}.pi_x = sqrt(pi_x^2 + pi_x_r^2);
        error{3}.pi_t = sqrt(pi_t^2 + pi_t_r^2);
        error{3}.pi_y = sqrt(pi_y^2 + pi_y_r^2);
        error{3}.pi_xty = sqrt(pi_x^2 + pi_x_r^2 + pi_t^2 + pi_t_r^2 + pi_y^2 + pi_y_r^2);
    else
        pi_x_r = nan;
        pi_t_r = nan;
        pi_y_r = nan;
        error{3}.pi_x_r = nan;
        error{3}.pi_t_r = nan;
        error{3}.pi_y_r = nan;
        error{3}.pi_y_mi_r = nan;
        error{3}.y_mi_r = nan;
        error{3}.eta_x_z_T = nan;
        error{3}.eta_u_z = nan;
        error{3}.eta_t_z_du = nan;
        error{3}.eta_t_z_u_w = nan;
        error{3}.eta_t_z_eta_w = nan;

        error{3}.pi_x = nan;
        error{3}.pi_t = nan;
        error{3}.pi_y = nan;
        error{3}.pi_xty = nan;
    end
    fprintf('done (%f seconds)\n', toc())

    %% ADAPTIVITY
    E = set_threshold_E(error_bk, params);
    %% REJECT
    if error_bk.pi_I > E

        %% Parametric
        fprintf('Computing reduced margin error indicators...'); tic();
        [pi_I_alpha, work, RM] = compute_pi_I_alpha(I_star, I, params, precompute);
        error_bk.pi_I_alpha = pi_I_alpha;
        error{2} = error_bk;
        J = mark_I(pi_I_alpha, work, RM, E, params);
        fprintf('marked indices...')
        disp(J);
        data_table = post_process(data_table, t + delta_t, delta_t, Z_I_star,  Z_I_star_lofi, ...
            I_star, I, E, error, ...
            work, RM, J, problem, params, fem, precompute, ...
            axArray);
        fprintf('done (%f seconds)\n', toc())

        fprintf('Refine approximation structure...'); tic();
        [Z_I_star, I_new, I_star_new] = refine_approximation(Z_I_star, I, I_star, J, problem, params, fem, 1);
        [Z_I_star_lofi, ~, ~] = refine_approximation(Z_I_star_lofi, I, I_star, J, problem, params, fem, 0);

        I = I_new;
        I_star = I_star_new;
        fprintf('done (%f seconds)\n', toc())

        fprintf('Recompute L2 norms...'); tic();
        precompute = compute_lagrange_norms(I_star, I, precompute);
        fprintf('done (%f seconds)\n', toc())

        delta_t = delta_t/params.k_shrink;
        %% Accept
    else
        fprintf('Accept step forward to time t=%f\n-----------------------\n',t+delta_t)
        pi_I_alpha = []; work = []; RM = []; J = [];

        pi_x = sqrt(pi_x^2 + pi_x_r^2);
        pi_t = sqrt(pi_t^2 + pi_t_r^2);
        pi_y = sqrt(pi_y^2 + pi_y_r^2);

        data_table = post_process(data_table, t + delta_t, delta_t, Z_I_star, Z_I_star_lofi,...
            I_star, I, E, error,...
            work, RM, J, problem, params, fem, precompute,...
            axArray);

        reference = compute_on_reference(t, delta_t, data_table, reference, Z_I_star, I_star, I, precompute, params, fem);
        
        t = t + delta_t;
        delta_t = delta_t*params.k_grow;
    end
    data_table(end,:)
end


if params.reference == 1
    reference.fem = fem;
end
%     save('reference.mat',"reference");

% end
