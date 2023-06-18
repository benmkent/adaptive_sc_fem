%% Generate path
addpath('./subroutines')
addpath('./tools')
addpath('./ifiss_modified')
addpath('./sparse-grid')
addpath('./bk')

addpath(genpath('../ifiss3.6'))
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
%     precompute = precompute_uQv(Z_I_star, precompute);
    precompute = precompute_interp_and_l2spatial(Z_I_star, I_star, I, fem, precompute);

    if params.jomp_estimator == 1
        fprintf('Compute BK hierarchical estimators...\n')
        [error_bk, Z_I_star, precompute] = error_estimate_bk(Z_I_star, Z_I_star_lofi, I_star, I, params, fem, precompute);
    else
        error_bk.pi = nan;
        error_bk.pi_I = nan;
        error_bk.pi_delta = inf;
        error_bk.pi_I_delta= inf;
        error_bk.pi_I_alpha = {};
    end
    error{2} = error_bk;

    fprintf('done (%f seconds)\n', toc())

    %% ADAPTIVITY
    E = set_threshold_E(error_bk, params);
    %% REJECT
    if error_bk.pi_I > E
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
