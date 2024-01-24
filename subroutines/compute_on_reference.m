function reference = compute_on_reference(t, delta_t, reference, Z_I_star, I_star, I, precompute, params, fem)
%COMPUTE_ON_REFERENCE Computes the approximation at the reference times and
%computes errors wrt reference approximation
%
% Inputs    t          sync time
%           delta_t     algorithm timestep
%           reference   reference data structure
%           Z_I_star    Collocation point data structure
%           I_star      Enhanced sparse grid
%           I           Sparse grid
%           precompute  precomputed data
%           params      approximation parameters
%           fem         finite element approximation options
%
% Outputs   reference   updated reference data structure


%% Set up
I_star_r = reduce_sparse_grid(I_star);
I_r = reduce_sparse_grid(I);

% If saving as reference then enhanced grid is used, else use standard grid
if params.reference == 1
    Iref = I_star;
    Iref_r = I_star_r;
    inds_z = 1:I_star_r.size;
else
    Iref = I;
    Iref_r = I_r;
    if I_star_r.size == I_r.size
        inds_z = 1:I_star_r.size
    else
        [~, inds_z, ~, ~] = compare_sparse_grids(I_star, I_star_r, I, I_r);
    end
end

% Find reference times in current interval
reference_times = reference.times;
inds = find((reference_times >= t) & (reference_times < t + delta_t));
reference_times = reference_times(inds);

% Collect solution at reference times
for ii = 1:Iref_r.size
    u_z = Z_I_star{inds_z(ii)}.u_z;
    t_z = Z_I_star{inds_z(ii)}.t_z;
    u_z_ref(ii,:,1:length(reference_times)) = interp1([t_z], [u_z]', reference_times)';
end

% Calculate expected value and approximate standard deviation
exp_u = zeros(size(u_z_ref,2), length(reference_times));
std_u=exp_u;
for ii = 1:length(reference_times)
    u_z_interpolated = interpolate_on_sparse_grid(Iref,Iref_r,[u_z_ref(:,:,ii)]', I_star_r.knots);

    exp_u(:,ii) = sum([u_z_ref(:,:,ii)]' .* Iref_r.weights,2);
    std_u(:,ii) = abs(sqrt(squeeze(sum([u_z_interpolated(:,:)].^2 .* I_star_r.weights,2)) - exp_u(:,ii).^2));
end

%% Compute error wrt reference approximation
if isempty(reference.reference_table)
    % If there is no reference approximation save errors as NaN
    ref_error(1:length(reference_times)) = nan;
    ref_error_H(1:length(reference_times)) = nan;
    global_ts_error(1:length(reference_times),Iref_r.size) = nan;
else
    % Reference error exists
    reference_table = reference.reference_table;
    notbound = Z_I_star{1}.notbound;

    % Compute L2 error and H1 error
    Q = fem.mass;
    H = fem.stiffness;
    fem_ref = reference.fem;
    Qref = fem_ref.mass;
    Href = fem_ref.stiffness;

    % For each reference time in current interval
    for ii_t = 1:length(reference_times)
        % Set up
        ref_row = reference_table(inds(ii_t),:);
        ref_I = ref_row.I{1};

        U = squeeze(u_z_ref(:,:,ii_t))';
        ref_U = ref_row.u{1};
        ref_I_r = reduce_sparse_grid(ref_I);

        u1Qv1 = zeros(size(ref_U,2));
        u2Qv2 = zeros(I_r.size);
        u1Qv2 = zeros(size(ref_U,2),I_r.size);
        u1Hv1 = zeros(size(ref_U,2));
        u2Hv2 = zeros(I_r.size);
        u1Hv2 = zeros(size(ref_U,2),I_r.size);

        % Remesh approximation onto (same or finer) reference mesh
        for ii_remesh = 1 : size(U,2)
            if strcmp(params.grid,'p1')
                if params.reference == 1
                    U_remesh = U;
                else
                    U_remesh(:,ii_remesh) = scattered_interpolant_bk(fem.T, fem.xy, U(:,ii_remesh), fem_ref.xy(:,1),fem_ref.xy(:,2));
                end
            elseif strcmp(params.grid,'q1')
                U_remesh(:,ii_remesh) = U(:,ii_remesh); %interp2(fem.xy(:,1),fem.xy(:,2), U(:,ii_remesh), fem_ref.xy(:,1),fem_ref.xy(:,2));
            else
                error();
            end
        end

        % Precompute weighted inner products on reference mesh
        for ii_Q = 1:size(ref_U,2)
            for ii_Q2 = 1:size(ref_U,2)
                u1Qv1(ii_Q,ii_Q2) = ref_U(:,ii_Q)' * Qref * ref_U(:,ii_Q2);
                u1Hv1(ii_Q,ii_Q2) = ref_U(:,ii_Q)' * Href * ref_U(:,ii_Q2);
                if ii_Q2 <= I_r.size
                    u1Qv2(ii_Q,ii_Q2) = ref_U(:,ii_Q)' * Qref * U_remesh(:,ii_Q2);
                    u1Hv2(ii_Q,ii_Q2) = ref_U(:,ii_Q)' * Href * U_remesh(:,ii_Q2);
                    if ii_Q <= I_r.size
                        u2Qv2(ii_Q,ii_Q2) = U(:,ii_Q)' * Q * U(:,ii_Q2);
                        u2Hv2(ii_Q,ii_Q2) = U(:,ii_Q)' * H * U(:,ii_Q2);
                    end
                end
            end
        end

        % Compute L^2_\rho(\Gamma) norms with precomputed data
        ref_I_r = sparse_grid_map_to_one_d_polynomials(ref_I,ref_I_r);
        I_r = sparse_grid_map_to_one_d_polynomials(I,I_r);
        ref_error(ii_t) = calc_l2rho_diff_quick_precalc(ref_I_r, I_r, u1Qv1, u2Qv2, u1Qv2, precompute.lagrange_product_integrals);
        ref_error_H(ii_t) = calc_l2rho_diff_quick_precalc(ref_I_r, I_r, u1Hv1, u2Hv2, u1Hv2, precompute.lagrange_product_integrals);

        %% Estimate global timestepping error at each collocation point
        % Find points I in ref_I
        if ref_I_r.size == I_r.size
            I_in_ref_I = 1:I_r.size;
            I_in_I = 1:I_r.size;
        else
            [~,I_in_ref_I,I_in_I,~] = compare_sparse_grids(ref_I, ref_I_r, I, I_r);
        end
        for ii_k = 1:length(I_in_ref_I)
            e_k = ref_U(:,I_in_ref_I(ii_k)) - U_remesh(:,I_in_I(ii_k));
            global_ts_error(ii_t,ii_k) = sqrt(e_k.' * Qref * e_k);
        end
    end
end

%% Save data to reference table
% Prepare data for each reference time in interval
ref_data = [];
for ii = 1:length(reference_times)
    row = table(reference_times(ii), {squeeze(u_z_ref(:,:,ii))'}, {Iref}, ref_error(ii), ref_error_H(ii), {global_ts_error(ii,:)},{exp_u(:,ii)}, {std_u(:,ii)},...
        'VariableNames',{'t','u','I','error','error_H','ts_error','exp','std'});
    ref_data = [ref_data; row];
end
% Append new reference data to reference.data_table
reference.data_table = [reference.data_table; ref_data];