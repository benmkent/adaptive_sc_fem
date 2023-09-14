function reference = compute_on_reference(t, delta_t, data_table, reference, Z_I_star, I_star, I, precompute, params, fem)

I_star_r = reduce_sparse_grid(I_star);
I_r = reduce_sparse_grid(I);

if params.reference == 1
    Iref = I_star;
    Iref_r = I_star_r;
    inds_z = 1:I_star_r.size;
else
    Iref = I;
    Iref_r = I_r;
    [~, inds_z, ~, ~] = compare_sparse_grids(I_star, I_star_r, I, I_r);
end

% Find reference times in interval
reference_times = reference.times;
inds = find((reference_times >= t) & (reference_times < t + delta_t));
reference_times = reference_times(inds);

% Compile solution at reference times
for ii = 1:Iref_r.size
    u_z_tplusdt = Z_I_star{inds_z(ii)}.u_z_tplusdt;
    t_z_tplusdt = Z_I_star{inds_z(ii)}.tplusdt;
    u_z = Z_I_star{inds_z(ii)}.u_z;
    t_z = Z_I_star{inds_z(ii)}.t_z;
%     u_z_ref(ii,:,1:length(reference_times)) = interp1([t_z_tplusdt, t_z], [u_z_tplusdt, u_z]', reference_times)';
    u_z_ref(ii,:,1:length(reference_times)) = interp1([t_z], [u_z]', reference_times)';
end

% Calculate exp and std
for ii = 1:length(reference_times)
    u_z_interpolated = interpolate_on_sparse_grid(Iref,Iref_r,[u_z_ref(:,:,ii)]', I_star_r.knots);

    exp_u(:,ii) = sum([u_z_ref(:,:,ii)]' .* Iref_r.weights,2);
    std_u(:,ii) = abs(sqrt(squeeze(sum([u_z_interpolated(:,:)].^2 .* I_star_r.weights,2)) - exp_u(:,ii).^2));
end

if isempty(reference.reference_table)
    ref_error(1:length(reference_times)) = nan;
    ref_error_H(1:length(reference_times)) = nan;
    global_ts_error(1:length(reference_times),Iref_r.size) = nan;
else
    reference_table = reference.reference_table;
    notbound = Z_I_star{1}.notbound;
    
    %% Compute L2 error and H1 error
    Q = fem.mass;
    H = fem.stiffness;
    fem_ref = reference.fem;
    Qref = fem_ref.mass;
    Href = fem_ref.stiffness;

    for ii_t = 1:length(reference_times)
        ref_row = reference_table(inds(ii_t),:);
        ref_I = ref_row.I{1};

        U = squeeze(u_z_ref(:,:,ii_t))';
        ref_U = ref_row.u{1};
        ref_I_r = reduce_sparse_grid(ref_I);

        %         u1Qv1 = uQv_ref{ii_t};

        u1Qv1 = zeros(size(ref_U,2));
        u2Qv2 = zeros(I_r.size);
        u1Qv2 = zeros(size(ref_U,2),I_r.size);
        u1Hv1 = zeros(size(ref_U,2));
        u2Hv2 = zeros(I_r.size);
        u1Hv2 = zeros(size(ref_U,2),I_r.size);

        for ii_remesh = 1 : size(U,2)
%             F1 = scatteredInterpolant(fem.xy(:,1),fem.xy(:,2),U(:,ii_remesh));
%             U_remesh(:,ii_remesh) = F1(fem_ref.xy(:,1),fem_ref.xy(:,2));
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



        for ii_Q = 1:size(ref_U,2)
            for ii_Q2 = 1:size(ref_U,2)            
                u1Qv1(ii_Q,ii_Q2) = ref_U(:,ii_Q)' * Qref * ref_U(:,ii_Q2);
                u1Hv1(ii_Q,ii_Q2) = ref_U(:,ii_Q)' * Href * ref_U(:,ii_Q2);
                if ii_Q2 <= I_r.size
                    u1Qv2(ii_Q,ii_Q2) = ref_U(:,ii_Q)' * Qref * U_remesh(:,ii_Q2);
                    u1Hv2(ii_Q,ii_Q2) = ref_U(:,ii_Q)' * Href * U_remesh(:,ii_Q2);
                    if ii_Q <= I_r.size
                        u2Qv2(ii_Q,ii_Q2) = U(:,ii_Q)' * Q * U(:,ii_Q2);
%                         u2Qv2_remesh(ii_Q,ii_Q2) = U_remesh(:,ii_Q)' * Qref * U_remesh(:,ii_Q2);
                        u2Hv2(ii_Q,ii_Q2) = U(:,ii_Q)' * H * U(:,ii_Q2);
%                         u2Hv2_remesh(ii_Q,ii_Q2) = U_remesh(:,ii_Q)' * Href * U_remesh(:,ii_Q2);
                    end
                end
            end
        end
        
%         u2Qv2 = u2Qv2_remesh;
%         u2Hv2 = u2Hv2_remesh;

        ref_I_r = sparse_grid_map_to_one_d_polynomials(ref_I,ref_I_r);
        I_r = sparse_grid_map_to_one_d_polynomials(I,I_r);
        ref_error(ii_t) = calc_l2rho_diff_quick_precalc(ref_I_r, I_r, u1Qv1, u2Qv2, u1Qv2, precompute.lagrange_product_integrals);
        ref_error_H(ii_t) = calc_l2rho_diff_quick_precalc(ref_I_r, I_r, u1Hv1, u2Hv2, u1Hv2, precompute.lagrange_product_integrals);

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

ref_data = [];
for ii = 1:length(reference_times)
    row = table(reference_times(ii), {squeeze(u_z_ref(:,:,ii))'}, {Iref}, ref_error(ii), ref_error_H(ii), {global_ts_error(ii,:)},{exp_u(:,ii)}, {std_u(:,ii)},...
        'VariableNames',{'t','u','I','error','error_H','ts_error','exp','std'});
    ref_data = [ref_data; row];
end

reference.data_table = [reference.data_table; ref_data];

% %% Plot
% figure(4); cla(); hold on;
% plot(data_table,'t','pi');
% plot(reference.data_table,'t','error');
% legend('show'); set(gca,'XScale','log','YScale','log');
% 
% figure(5); hold on;
% for ii = 1:size(ref_data,1)
%     t_plot = ref_data{ii,'t'};
%     ts_errors = ref_data{ii,'ts_error'};
%     scatter(t_plot,ts_errors{1},'k.');
% end
% t_plot = Z_I_star(1).tplusdt;
% ts_errors = [Z_I_star(:).ge_estimate];
% scatter(t_plot,ts_errors,'b.');
% set(gca,'XScale','log','YScale','log');
% drawnow();
% end