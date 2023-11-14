function [data_table,fem,problem,params,reference] = adaptive_sc_fem(problem, params, reference)
%ADAPTIVE_SC_FEM Construct an adaptive in time SC-FEM approximation
% Uses Algorithm 1 from thesis.
% Inputs
%   problem     A structure detailing the advection--diffusion problem
%   params      A structure detailing the approximation parameters
%   reference   A structure containing a reference approximation and times

%% Initialise
% Line 3: Construct the MI sets $I$ and $I^*$.
fprintf('Initialise sparse grid...'); tic();
[I_star, I] = initialise_sparse_grid(problem, params);
I_star_r = reduce_sparse_grid(I_star);
fprintf('done (%f seconds)\n', toc())

% Line 4: Construct FEM matrices
fprintf('Initialise FEM matrices...'); tic();
% First construct the mesh $M^h$
fem = initialise_fem_mesh(params);
% Then construct the FEM matrices $D^h(z)$, $W^h(z)$, $f^h(z)$.
fem = initialise_fem_matrices(problem,params,fem);
fprintf('done (%f seconds)\n', toc());
fprintf('Total DOF %d\n', size(fem.xy,1) * I_star_r.size);

% Initialise a structure to contain the collocation points and associated
% data.
% Line 6 included.
fprintf('Initialise approximation structure...'); tic();
Z_I_star = initialise_approximation_structure(params, fem, I_star);
% For hierarchical error estimator we require a lofi approximation for
% estimating the global timestepping error.
if strcmp(params.adapt_type,'hierarchical')
    params_lofi = params;
    params_lofi.letol = params.letol_lofi;
    params_lofi.dt0 = params.dt0_lofi;

    Z_I_star_lofi = initialise_approximation_structure(params_lofi, fem, I_star);
end
fprintf('done (%f seconds)\n', toc())

% Initialise a structure for precompute data
fprintf('Initialise precomputed data (Lagrange product integrals)...'); tic();
precompute = initialise_precomputed_data(I_star, I, params);
fprintf('done (%f seconds)\n', toc())

delta_t = params.t0;
t_r=0;

%% Initialise datatable and figures for storing data
data_table = write_to_data_table(params);
axArray = initialise_figures();

if strcmp(params.adapt_type,'residual')
    error_residual.pi_x = 0;
    error_residual.pi_t = 0;
    error_residual.pi_y = 0;
    error_residual.pi_xty = 0;
end
%% Start adaptive loop
% Line 8
while t_r < params.T
    %% SOLVE
    % Line 9: for each collocation point propagate forward to $t+\tau$.
    fprintf('Propagate SC point solutions to at least t=%f...\n',t_r+delta_t); tic();
    for ii_z = 1:length(Z_I_star)
        fprintf('.... (target t= %f) %d of %d',t_r+delta_t, ii_z, length(Z_I_star));
        % Line 10: Adaptive timestepping
        % Line 11,12: Update data in structure
        Z_I_star{ii_z} = propagate_sc_point( Z_I_star{ii_z}, t_r + delta_t, params);
        if strcmp(params.adapt_type,'hierarchical')
            % Compute approximations with LoFi timestepping
            Z_I_star_lofi{ii_z} = propagate_sc_point(Z_I_star_lofi{ii_z},...
                t_r + delta_t, params_lofi);
        end
        dt_tplusdt = Z_I_star{ii_z}.dt_z_tplusdt;
        fprintf('...(%f timestep at target t)', dt_tplusdt);
        fprintf('...done\n');
    end
    fprintf('done (%f seconds)\n', toc())

    %% ESTIMATE
    fprintf('Compute error estimates...'); tic();
    % Precompute the the $u_i' Q v_j$ $u_i' A v_j$ products to speed up computing norms of
    % the form $\Vert \sum_{i} u_i L_i - \sum_{j} u_j L_j \Vert$
    fprintf('.... Precompute vector-matrix-vector products\n')
    precompute = precompute_interp_and_l2spatial(Z_I_star, I_star, I, fem, precompute);

    % Line 15: Compute the error estimators
    switch params.adapt_type
        case 'hierarchical'
            % See Chapter 3 or Efficient Adaptive Stochastic Collocation
            % Strategies for Advection–Diffusion Problems with Uncertain
            % Inputs
            fprintf('Compute BK hierarchical estimators...\n')
            [error_bk, Z_I_star] = error_estimate_bk(Z_I_star, Z_I_star_lofi, ...
                I_star, I, params, precompute);
            error{2} = error_bk;
        case 'residual'
            % See Chapter 5
            fprintf('Compute residual estimators\n')
            [error_residual] = error_estimate_residual(error_residual, t_r, delta_t,...
                Z_I_star, I_star, I, params, fem, precompute);
    end
    fprintf('done (%f seconds)\n', toc())

    %% ADAPT
    switch params.adapt_type
        case 'hierarchical'
            % Line 16: Set tolerance $E = c_tol * \pi^{I^*}$
            E = set_threshold_E(error_bk, params);
            % Line 18: Test $E > \pi^{I}$.
            adapt_flag = error_bk.pi_I > E;
        case 'residual'
            % Set threshold for residual based algorithm
            E = max([error_residual.pi_t_r, 1e0]);
%             E = 10*delta_t;
            % Test against threshold
            adapt_spatial = error_residual.pi_x_r > E;
            adapt_param = error_residual.pi_y_r > E;
            adapt_flag =  adapt_spatial || adapt_param;
        otherwise
            % If no adaptivity set tolerance to $\infty$
            E = inf;
    end

    %% REJECT
    if adapt_flag == 1
        switch params.adapt_type
            case 'hierarchical'
                %% Parametric
                % Line 19: Compute the error indicators in 19.
                fprintf('Computing reduced margin error indicators...'); tic();
                [pi_I_alpha, work, RM] = compute_pi_I_alpha(I_star, I, params, precompute);
                error_bk.pi_I_alpha = pi_I_alpha;
                % Line 20: Mark subset of multi-indices
                J = mark_I(pi_I_alpha, RM, params);
                fprintf('marked indices...')
                disp(J);
                % Write data to data_table
                data_table = post_process(data_table, t_r + delta_t, delta_t, Z_I_star,  Z_I_star_lofi, ...
                    I_star, I, E, error_bk,...
                    work, RM, J, problem, params, fem, precompute, ...
                    axArray);
                fprintf('done (%f seconds)\n', toc())

                % Line 21-24: Initialise approximation structure
                fprintf('Refine approximation structure...'); tic();
                [Z_I_star, I_new, I_star_new] = refine_approximation(tr, Z_I_star, I, I_star, J, problem, params, fem, 1);
                [Z_I_star_lofi, ~, ~] = refine_approximation(tr, Z_I_star_lofi, I, I_star, J, problem, params, fem, 0);

                % Line 25-26: Update multi-index sets
                I = I_new;
                I_star = I_star_new;
                % Line 27: Shrink algorithm timestep.
                delta_t = delta_t/params.k_shrink;
                fprintf('done (%f seconds)\n', toc())

                % Recompute $L_{\rho}^2(\Gamma)$ norms for new sets of
                % interpolation polynomials.
                fprintf('Recompute L2 norms...'); tic();
                precompute = compute_lagrange_norms(I_star, I, precompute);
                fprintf('done (%f seconds)\n', toc())
            case 'residual'
                if adapt_spatial == 1
                    %% Spatial refinement
                    I_r = reduce_sparse_grid(I);
                    [elerr] = flatten_spatial_estimator(error_residual.eta_x_z_T,I_r);

                    plot_data_tifiss(1,1,Z_I_star{1}.u_z(:,1),elerr,fem.ev,fem.xy);

                    markstrat = 2; % dorfler
                    threshold = min([0.3,max([0.1,1.5*(error_residual.pi_x_r - E)/error_residual.pi_x_r])]);

                    % Mark subset of elements
                    [Mset] = marking_strategy_fa(elerr,markstrat,threshold);

                    % Refine FEM mesh for every collocation point
                    [Z_I_star,fem] = refine_fem_mesh(Z_I_star,Mset,fem,params);
                    fem = initialise_fem_matrices(problem,params,fem);
                    % Update data structure
                    [Z_I_star_refined, ~ ,~] = refine_approximation(t_r, Z_I_star, I, I_star, [], problem, params, fem, 1);
                    Z_I_star = Z_I_star_refined;
                end
                if adapt_param == 1
                    %% Parametric refinement
                    J = mark_I(error_residual.pi_y_mi_r, error_residual.eta_mi, params);
                    fprintf('marked indices...')
                    disp(J);
                    fprintf('done (%f seconds)\n', toc())

                    fprintf('Refine approximation structure...'); tic();
                    [Z_I_star, I_new, I_star_new] = refine_approximation(Z_I_star, I, I_star, J, problem, params, fem, 1);

                    I = I_new;
                    I_star = I_star_new;
                    fprintf('done (%f seconds)\n', toc())

                    % Only need L1 norms so could seperate out the
                    % following function.
                    fprintf('Recompute L2 norms...'); tic();
                    precompute = compute_lagrange_norms(I_star, I, precompute);
                    fprintf('done (%f seconds)\n', toc())
                else
                    J=[];
                end
                % Write to data table
                data_table = post_process(data_table, t_r + delta_t, delta_t, Z_I_star,  [], ...
                    I_star, I, E, error_residual, ...
                    [], error_residual.y_mi_r, J, problem, params, fem, precompute, ...
                    axArray);
        end
        %% Accept
    else
        fprintf('Accept step forward to time t=%f\n-----------------------\n',t_r+delta_t)
        switch params.adapt_type
            case 'hierarchical'
                work = []; RM = []; J = [];
                % Write data to data_table
                data_table = post_process(data_table, t_r + delta_t, delta_t, Z_I_star, Z_I_star_lofi,...
                    I_star, I, E, error_bk,...
                    work, RM, J, problem, params, fem, precompute,...
                    axArray);
            case 'residual'
                % Accumulate total error
                work = []; RM = []; J = [];
                error_residual.pi_x = sqrt(error_residual.pi_x^2 + error_residual.pi_x_r^2);
                error_residual.pi_t = sqrt(error_residual.pi_t^2 + error_residual.pi_t_r^2);
                error_residual.pi_y = sqrt(error_residual.pi_y^2 + error_residual.pi_y_r^2);

                % Write data to data_table
                data_table = post_process(data_table, t_r + delta_t, delta_t, Z_I_star, [],...
                    I_star, I, E, error_residual,...
                    work, RM, J, problem, params, fem, precompute,...
                    axArray);
        end


        % Comptue error with respect to precomputed reference approx and
        % save out the approximation at the reference times.
        reference = compute_on_reference(t_r, delta_t, reference, Z_I_star, I_star, I, precompute, params, fem);

        % Line 29/32: Update the algorithm timestep and time t.
        t_r = t_r + delta_t;
        delta_t = delta_t*params.k_grow;
    end
    %data_table(end,:)
end

if params.reference == 1
    reference.fem = fem;
end
end
