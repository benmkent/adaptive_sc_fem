function [eta_y,eta_y_mi,eta_mi] = compute_interpolation_errors(tr,trp1,Z_I_star, I, I_star, fem, params, precompute)
%COMPUTE_INTERPOLATION_ERRORS Compute thte interpolation error estimators
%
% Inputs    tr          time interval start
%           trp1        time interval end
%           Z_I_star    data structure of collocation points
%           I           Sparse grid
%           I_star      Enhanced sparse grid
%           fem         Finite element structure
%           params      Approximation parameters
%           precompute  Precomputed data
%
% Outputs   eta_y       interpolation error estimator
%           eta_y_mi    interpolation error estimator for each MI
%           eta_mi      corresponding MI for each estimate

%% Set up params
ev = fem.ev;
xy = fem.xy;
x = xy(:,1);
y= xy(:,2);
I_star_r = reduce_sparse_grid(I_star);
I_r = reduce_sparse_grid(I);
if I_star_r.size==I_r.size
    pts_in_both_star = 1:I_star_r.size;
else
[pts_in_only_star, pts_in_both_star, pts_in_both, pts_in_only ] = compare_sparse_grids(I_star,I_star_r,I,I_r);
end

Z_I = Z_I_star(pts_in_both_star);

%% Precompute shared timesteps
% Union of the timegrids for each collocation point.
% The approximation is consequently piecewise linear on this.
% Find a common set of times for integrating.
t_common=[];
for ii = 1:length(pts_in_both_star)
    t_common = [t_common; Z_I_star{pts_in_both_star(ii)}.t_z];
end
t_common = unique([tr;trp1;t_common]); % finds unique AND sorts
% Restrict to times in [t_r=t,...,t_r+1=t+delta_t
t_common(t_common < tr) = [];
t_common(t_common > trp1) = [];

% Interpolate each solution to common times
for ii_z = 1:length(Z_I)
    Z_I{ii_z}.u_t_common = interp1(Z_I{ii_z}.t_z, (Z_I{ii_z}.u_z).', t_common).';
end
dt = diff(t_common);
n_k = length(t_common);

%% Now compute estimators
% Now define suitable batches for estimating interpolation error
% See Remark 6.1 Guignard and Nobile 2018
% We sequentially form interpolants by considering only considering the RM
% (ensuring admissiblity of ((G \cup RM) \cup RM_2) \cup RM_3)...
% where the original margin is M = \CUP RM_i

% Identify margin
layer_rm = 1;
miset = get_mi_set(I);
G{layer_rm} = miset;
[RM{layer_rm}, M] = sparse_grid_reduced_margin(G{layer_rm});
nM = size(M,1);
M_remaining = setdiff(M,  [G{layer_rm};RM{layer_rm}]);

% Identify MI in margin that are in RM of G at each "layer"
while ~isempty(M_remaining)
    G{layer_rm+1} = [G{layer_rm};RM{layer_rm}];
    layer_rm = layer_rm +1;
    % Identify Rm for current layer
    [RM{layer_rm}, ~] = sparse_grid_reduced_margin(G{layer_rm});

    % Only keep RM that is in original margin
    RM{layer_rm} = intersect(RM{layer_rm}, M,'rows');

    % Identify remaining indices
    M_remaining = setdiff(M, [G{layer_rm};RM{layer_rm}]);
end

%% Compute advection and diffusion fields for all collocation points
% Advection and diffusion fields are evalauted at a set of quadrature
% points on each element for each collocation point

% Define a superset of all points for interpolation
[~,C] = check_set_admissibility([miset;M]);
[Isuper] = smolyak_grid_multiidx_set(C,params.knot_fn,params.lev2knots,I);
Isuper_r = reduce_sparse_grid(Isuper);

% Define quadature points in space
% Construct 2D gaussian rule over the reference triangle
nngpt = 7;
[s,t,wt] = triangular_gausspoints(nngpt);
% For each gauss pt define the x and y pts in each element
for ivtx = 1:3
    xl_v(:,ivtx) = x(ev(:,ivtx));
    yl_v(:,ivtx) = y(ev(:,ivtx));
end
nel = size(xl_v,1);
for igpt = 1:nngpt
    sigpt=s(igpt);
    tigpt=t(igpt);
    xl_v(:,ivtx) = x(ev(:,ivtx));
    yl_v(:,ivtx) = y(ev(:,ivtx));
    [phi_e,dphids,dphidt] = tshape(sigpt,tigpt);
    xx = zeros(nel,1);
    yy=xx;
    for ivtx=1:3
        xx = xx + phi_e(ivtx) * xl_v(:,ivtx);
        yy = yy + phi_e(ivtx) * yl_v(:,ivtx);
    end
    [jac_igpt(:,igpt), invjac_gpt(:,igpt),~,~,~] = tderiv(sigpt,tigpt,xl_v,yl_v);

    x_gpt(:,igpt) = xx;
    y_gpt(:,igpt) = yy;

    % Now evaluate diff and wind field at quadature points for z in Z^I
    for ii_z = 1:Isuper_r.size;
        diff_on_super(:,ii_z,igpt) = fem.a_fn(xx,yy,Isuper_r.knots(:,ii_z))*[1;0;0;0];
        wind_evaluate = fem.wind_fn(xx,yy,Isuper_r.knots(:,ii_z));
        windx_on_super(:,ii_z,igpt) = wind_evaluate(:,1);
        windy_on_super(:,ii_z,igpt) = wind_evaluate(:,2);
    end
end

%% Compute the interpolation error indicators at each discrete time
% Now for each time interval
for ii_k = 1:(n_k)
    fprintf('... for common timestep %d of %d\n',ii_k,n_k)
    % Compute spatial gradients for this timestep
    U_I = zeros(size(fem.xy,1),length(Z_I));
    dUx = zeros(size(fem.ev,1),length(Z_I));
    dUy = zeros(size(fem.ev,1),length(Z_I));
    for ii_z = 1:length(Z_I)
        U_I(:,ii_z) = [Z_I{ii_z}.u_t_common(:,ii_k)];
        [dUx(:,ii_z), dUy(:,ii_z)] = get_element_gradients(U_I(:,ii_z), fem);
    end

    % Interpolate gradients to superset of collocation points
    dUx_super = interpolate_on_sparse_grid(I,I_r, dUx, Isuper_r.knots);
    dUy_super = interpolate_on_sparse_grid(I,I_r, dUy, Isuper_r.knots);

    % Multiply by diff and advection fields
    % Note this assumes scalar valued diffusion!
    for igpt = 1:nngpt
        aDxU_super(:,:,igpt) = diff_on_super(:,:,igpt) .* dUx_super;
        aDyU_super(:,:,igpt) = diff_on_super(:,:,igpt) .* dUy_super;
        wxDxU_super(:,:,igpt) = windx_on_super(:,:,igpt) .* dUx_super;
        wyDyU_super(:,:,igpt) = windy_on_super(:,:,igpt) .* dUy_super;
    end

    jj_mi=0;
    %% For each layer G
    % For each layer G we compute the error indicators for each MI in the
    % corresponding RM
    for ii_l = 1:length(G)
        jj_mi = jj_mi+1;
        % Interpolate terms from grid I to grid G
        Gii = G{ii_l};
        I_G = smolyak_grid_multiidx_set(Gii,params.knot_fn,params.lev2knots,I);
        I_G_r = reduce_sparse_grid(I_G);

        [pts_in_super_only, pts_in_G_in_super] = compare_sparse_grids(Isuper,Isuper_r,I_G,I_G_r);
        aDxU_G = aDxU_super(:,pts_in_G_in_super,:);
        aDyU_G = aDyU_super(:,pts_in_G_in_super,:);
        wxDxU_G = wxDxU_super(:,pts_in_G_in_super,:);
        wyDyU_G = wyDyU_super(:,pts_in_G_in_super,:);

        % Interpolate from G to super points and subtract
        for igpt=1:nngpt
            aDxU_G_on_super(:,:,igpt) = interpolate_on_sparse_grid(I_G,I_G_r, aDxU_G(:,:,igpt), Isuper_r.knots);
            aDyU_G_on_super(:,:,igpt) = interpolate_on_sparse_grid(I_G,I_G_r, aDyU_G(:,:,igpt), Isuper_r.knots);
            wxDxU_G_on_super(:,:,igpt) = interpolate_on_sparse_grid(I_G,I_G_r, wxDxU_G(:,:,igpt), Isuper_r.knots);
            wyDyU_G_on_super(:,:,igpt) = interpolate_on_sparse_grid(I_G,I_G_r, wyDyU_G(:,:,igpt), Isuper_r.knots);

            aDxU_super_minus_G(:,:,igpt) = aDxU_super(:,:,igpt)  - aDxU_G_on_super(:,:,igpt);
            aDyU_super_minus_G(:,:,igpt) = aDyU_super(:,:,igpt)  - aDyU_G_on_super(:,:,igpt);
            wxDxU_super_minus_G(:,:,igpt) = wxDxU_super(:,:,igpt)  - wxDxU_G_on_super(:,:,igpt);
            wyDyU_super_minus_G(:,:,igpt) = wyDyU_super(:,:,igpt)  - wyDyU_G_on_super(:,:,igpt);
        end

        % Now we know for pts in G the differnece is zero
        % Compute L^2(D) norms for products of the differences
        % dUQdU(pts_in_both_new) = 0;
        % Compute for the pts in super grid only.
        l2pairs_diff_gpt = zeros(length(pts_in_super_only),length(pts_in_super_only),nel,nngpt);
        l2pairs_wind_gpt = zeros(length(pts_in_super_only),length(pts_in_super_only),nel,nngpt);
        for ii = 1:length(pts_in_super_only)
            for jj = 1:ii
                for igpt = 1:nngpt
                l2pairs_diff_gpt(ii,jj,:,igpt) = ...
                    (aDxU_super_minus_G(:,(pts_in_super_only(ii)),igpt) .* ...
                    aDxU_super_minus_G(:,(pts_in_super_only(jj)),igpt)) +...
                    (aDyU_super_minus_G(:,(pts_in_super_only(ii)),igpt) .* ...
                    aDyU_super_minus_G(:,(pts_in_super_only(jj)),igpt)).*...
                    jac_igpt(:,igpt) .* jac_igpt(:,igpt).* ...
                    wt(igpt);
                l2pairs_wind_gpt(ii,jj,:,igpt) = ((wxDxU_super_minus_G(:,(pts_in_super_only(ii)),igpt) + ...
                    wyDyU_super_minus_G(:,(pts_in_super_only(ii)),igpt)) .*...
                    (wxDxU_super_minus_G(:,(pts_in_super_only(jj)),igpt) + ...
                    wyDyU_super_minus_G(:,(pts_in_super_only(jj)),igpt))).*...
                    jac_igpt(:,igpt) .* jac_igpt(:,igpt).* ...
                    wt(igpt);
                end
            end
        end
        %% Sum over quadrature points on spatial domain
        % We sum over the elements and the gauss points
        l2pairs_diff = sum(l2pairs_diff_gpt,[3,4]);
        l2pairs_wind = sum(l2pairs_wind_gpt,[3,4]);
        % use symmetry to complete matrix
        l2pairs_diff = l2pairs_diff + l2pairs_diff.' - diag(diag(l2pairs_diff));
        l2pairs_wind = l2pairs_wind + l2pairs_wind.' - diag(diag(l2pairs_wind));
        
        % Now embed this in larger matrix for all collocation points
        l2pairs_diff_full = sparse(Isuper_r.size,Isuper_r.size);
        l2pairs_diff_full(pts_in_super_only,pts_in_super_only) = l2pairs_diff;
        l2pairs_wind_full = sparse(Isuper_r.size,Isuper_r.size);
        l2pairs_wind_full(pts_in_super_only,pts_in_super_only) = l2pairs_wind;        

        %% We are now in a position to integrate over parameter domain
        % For each multi-index in this layer
        RM_l = RM{ii_l};
        for ii_mi = 1:size(RM_l,1)
            mi = RM_l(ii_mi,:);
            fprintf('...... for  mi # %d of %d\n', jj_mi, nM);

            %% Build up sparse grid for G \cup mi
            Gmi = sortrows([Gii;mi],'ascend');
            I_G_mi = smolyak_grid_multiidx_set(Gmi,params.knot_fn,params.lev2knots);
            I_G_mi_r = reduce_sparse_grid(I_G_mi);
            [I_G_mi_r] = sparse_grid_map_to_one_d_polynomials(I_G_mi,I_G_mi_r);

            % Identify pts in super grid
            if I_G_mi_r.size ~= Isuper_r.size
                [ ~,inds_I_G_mi_in_Isuper,~,~] = compare_sparse_grids(Isuper,Isuper_r,I_G_mi,I_G_mi_r);
            else
                inds_I_G_mi_in_Isuper = 1:Isuper_r.size;
            end

            % Extract appropriate subset of precomputed l2 pairs 
            l2pairs_diff_mi = l2pairs_diff_full(inds_I_G_mi_in_Isuper,inds_I_G_mi_in_Isuper);
            l2pairs_wind_mi = l2pairs_wind_full(inds_I_G_mi_in_Isuper,inds_I_G_mi_in_Isuper);
            
            % Compute L^2_\rho(\Gamma) norms
            eta_diff(ii_k,jj_mi) = calc_l2rho_quick_precalc(I_G_mi_r,l2pairs_diff_mi, precompute.lagrange_product_integrals);
            eta_wind(ii_k,jj_mi) = calc_l2rho_quick_precalc(I_G_mi_r,l2pairs_wind_mi, precompute.lagrange_product_integrals);
            
            % Use known Poincare and compute estimates
            Poincare=fem.poincare;
            eta_y_k_mi(ii_k,jj_mi) = (eta_diff(ii_k,jj_mi) + Poincare*eta_wind(ii_k,jj_mi));
            eta_mi(jj_mi,:) = mi;
            jj_mi=jj_mi+1;
        end
    end
end

%% Integrate over time interval to give error indicators.
eta_y_mi = sqrt(jj_mi*2/3 * dt(:).' * (eta_y_k_mi(1:end-1,:).^2 + eta_y_k_mi(2:end,:).^2));
eta_y = sqrt(sum(eta_y_mi).^2);
