function fem = initialise_fem_matrices(problem, params)
%% Define geometry of spatial domain (using ifiss3.6)
pde=2; domain=1;
rsquare_domain_modified(params.grid_param, params.grid_type)
cd '../datafiles/'
load('square_grid.mat','grid_type','y','x','xy','mv','bound','mbound','outbc');
cd '../../adaptive_sc_fem'

[ev,ebound] = q1grid(xy,mv,bound,mbound);
switch problem.name
    case 'one-d'
        wind_fn = @(x,y,nel) wind_fn_parameterised(x,y,[0,0], 1);
        [A,N,Q,epe,eph,epw,n2] = femq1_cd_bk(xy,ev,wind_fn);
        N0 = N;


        epsilon = 1/10;
        sigma = 0.3;
        Diff = epsilon * A;
        % Conv = @(y) ((10^.5 + 10^-.5)/2 + y(1)*(10^.5-10^-.5)/2)*N0; % + y(:,2)*N2  + y(:,3)*N3 + y(:,4)*N4;
        Conv = @(y) N0 + sigma * y(1)*N0;

        AN = @(y) Diff + Conv(y); % Define the advection-diffusion matrix as a sum of fns xi of the parameters y
    case 'one-d-diff'
        wind_fn = @(x,y,nel) wind_fn_parameterised(x,y,[0,0], 1);
        [A,N,Q,epe,eph,epw,n2] = femq1_cd_modified(xy,ev,wind_fn);
        N0 = N;


        epsilon = 1/10;
        sigma = 0.3;
        Diff = @(y) (1+ sigma  * y(1)) * epsilon * A;
        % Conv = @(y) ((10^.5 + 10^-.5)/2 + y(1)*(10^.5-10^-.5)/2)*N0; % + y(:,2)*N2  + y(:,3)*N3 + y(:,4)*N4;
        Conv = @(y) N0;

        AN = @(y) Diff(y) + Conv(y); % Define the advection-diffusion matrix as a sum of fns xi of the parameters y

    case 'eigel'
        wind_fn = @(x,y,nel) wind_fn_parameterised(x,y,[0,0], 1);
        [A,N,Q,epe,eph,epw,n2] = femq1_cd_modified(xy,ev,wind_fn);
        N0 = N;
        A0 = A;

        epsilon = 0.1;
        nDim = problem.n;

        alpha_bar = 0.547;
        k = @(mm) floor(-0.5 + sqrt(1/4 + 2*mm));
        beta_1 = @(mm) mm - k(mm)*(k(mm)+1)/2;
        beta_2 = @(mm) k(mm) - beta_1(mm);
        alpha_m = @(mm) alpha_bar * mm^-2;

        xmap = @(x) 2*x - 1;
        for mm = 1:nDim
            diff_fn = @(x1,x2) alpha_m(mm) .* cos(2*pi* beta_1(mm) *x1) .* cos(2*pi*beta_2(mm)*x2);
            diff_fn_shifted = @(x1,x2) diff_fn(xmap(x1),xmap(x2));
            [Amm{mm},~,~,epe,eph,epw,n2] = femq1_cd_modified(xy,ev,wind_fn,diff_fn_shifted);
        end

        Diff = @(y) A0 + combine_matrices(y,Amm);
        % Conv = @(y) ((10^.5 + 10^-.5)/2 + y(1)*(10^.5-10^-.5)/2)*N0; % + y(:,2)*N2  + y(:,3)*N3 + y(:,4)*N4;
        Conv = @(y) N0;

        AN = @(y) epsilon*Diff(y) + Conv(y); % Define the advection-diffusion matrix as a sum of fns xi of the parameters y

    case 'quadrants-isotropic'
        %Wind fn to generate the matrix N
        sigma = [0.3,0.3,0.3,0.3];
        wind_fn = @(x,y,nel) wind_fn_parameterised(x,y,[0,0], 1);
        [A,N,Q,epe,eph,epw,n2] = femq1_cd_modified(xy,ev,wind_fn);
        N0 = N;
        wind_fn = @(x,y,nel) wind_fn_parameterised(x,y,[0.5,-0.5], 2);
        [A,N,Q,epe,eph,epw,n2] = femq1_cd_modified(xy,ev,wind_fn);
        N1 = N;
        wind_fn = @(x,y,nel) wind_fn_parameterised(x,y,[-0.5,-0.5], 2);
        [A,N,Q,epe,eph,epw,n2] = femq1_cd_modified(xy,ev,wind_fn);
        N2 = N;
        wind_fn = @(x,y,nel) wind_fn_parameterised(x,y,[-0.5,0.5], 2);
        [A,N,Q,epe,eph,epw,n2] = femq1_cd_modified(xy,ev,wind_fn);
        N3 = N;
        wind_fn = @(x,y,nel) wind_fn_parameterised(x,y,[0.5,0.5], 2);
        [A,N,Q,epe,eph,epw,n2] = femq1_cd_modified(xy,ev,wind_fn);
        N4 = N;

        Conv = @(y) N0 + sigma(1) * y(1)*N1 + sigma(2) * y(2)*N2 + sigma(3) * y(3)*N3 + sigma(4) * y(4)*N4;
        epsilon = 1/10;
        Diff = epsilon * A;
        AN = @(y) Diff + Conv(y); % Define the advection-diffusion matrix as a sum of fns xi of the parameters y
    case 'quadrants-anisotropic'
        %Wind fn to generate the matrix N
        sigma = [0.1,0.3,0.01,0.2];
        wind_fn = @(x,y,nel) wind_fn_parameterised(x,y,[0,0], 1);
        [A,N,Q,epe,eph,epw,n2] = femq1_cd_modified(xy,ev,wind_fn);
        N0 = N;
        wind_fn = @(x,y,nel) wind_fn_parameterised(x,y,[0.5,-0.5], 2);
        [A,N,Q,epe,eph,epw,n2] = femq1_cd_modified(xy,ev,wind_fn);
        N1 = N;
        wind_fn = @(x,y,nel) wind_fn_parameterised(x,y,[-0.5,-0.5], 2);
        [A,N,Q,epe,eph,epw,n2] = femq1_cd_modified(xy,ev,wind_fn);
        N2 = N;
        wind_fn = @(x,y,nel) wind_fn_parameterised(x,y,[-0.5,0.5], 2);
        [A,N,Q,epe,eph,epw,n2] = femq1_cd_modified(xy,ev,wind_fn);
        N3 = N;
        wind_fn = @(x,y,nel) wind_fn_parameterised(x,y,[0.5,0.5], 2);
        [A,N,Q,epe,eph,epw,n2] = femq1_cd_modified(xy,ev,wind_fn);
        N4 = N;

        Conv = @(y) N0 + sigma(1) * y(1)*N1 + sigma(2) * y(2)*N2 + sigma(3) * y(3)*N3 + sigma(4) * y(4)*N4;
        epsilon = 1/10;
        Diff = epsilon * A;
        AN = @(y) Diff + Conv(y); % Define the advection-diffusion matrix as a sum of fns xi of the parameters y

    case 'eddies'
        n_dim = problem.n;
        sigma = problem.sigma;
        sqrt_n = round(sqrt(n_dim));
        dist_centres = 2/sqrt_n;
        centres = -1 + dist_centres*0.5 + dist_centres * (0:(sqrt_n-1));

        wind_fn = @(x,y,nel) wind_fn_parameterised(x,y,[0,0], 1);
        [A,N,Q,epe,eph,epw,n2] = femq1_cd_modified(xy,ev,wind_fn);
        N0 = N;

        for ii= 1:sqrt_n
            for jj = 1:sqrt_n
                x_w = centres(ii);
                y_w = centres(jj);
                wind_fn = @(x,y,nel) wind_fn_parameterised(x,y,[x_w,y_w], sqrt_n);
                [A,N,Q,epe,eph,epw,n2] = femq1_cd_modified(xy,ev,wind_fn);
                Nij{(ii-1)*sqrt_n + jj} = sigma * N;
            end
        end
        Conv = @(y) N0 + combine_matrices(y, Nij); %N0 + sum(bsxfun(@times,reshape(y,1,1,[]),Nij),1);
        epsilon = 1/10;
        Diff = epsilon * A;
        AN = @(y) Diff + Conv(y); % Define the advection-diffusion matrix as a sum of fns xi of the parameters y

    case 'djs-div-free'
        wind_fn = @(x,y,nel) wind_fn_parameterised(x,y,[0,0], 1);
        [A,N0,Q,epe,eph,epw,n2] = femq1_cd_modified(xy,ev,wind_fn);

        %         Ben (Catherine, for information)
        % As discussed yesterday I am attaching a bunch of m-files
        % Note that these are designed to work with Q2 approximation not Q1
        %
        % ---------------------------------------------
        % Start with linstab_startdata.m lines 19 to 51
        % line 20 ...
        % xyp are the Q1 coordinates of the current grid defined by xy and mv
        % These are the vertices coordinates
        %[xyp,mp,map] = q2q1map(xy,mv);
        [ ~, xyp, boundp] = rsquare_domain_modified_fn(6,1);
        np = length(xyp);

        % lines 21 to 38
        % sets up discrete covariance matrix Cvar
        % (setting delta=1 generates rapidly decaying eigenvalues ...)

        % Set up covariance operator
        dx = abs( xyp(:,1)*ones(1,np) - ones(np,1)*xyp(:,1)' );
        dy = abs( xyp(:,2)*ones(1,np) - ones(np,1)*xyp(:,2)' );

        %corrL = 2.0;    %% reference length L for the domain
        %dist = sqrt( dx.^2 + dy.^2 );
        %delta = 1;
        %Cvar = exp(-(dist./corrL).^(1+delta));

        delta = 1; %problem.delta; % delta = 1;
        corrx = problem.corrlength; %2; % problem.corrx; %         corrx = 2;
        corry = problem.corrlength; %2; %problem.corry; %corry = 2;
        dist = sqrt((dx/corrx).^2 + (dy/corry).^2);
        Cvar = exp(-dist.^(1+delta));

        %dist = dx/corrx + dy/corry;
        %Cvar = exp(-dist);

        % line 43
        % cov is a multiplier ("coefficient of variation")

        cov   = 1e-1; %problem.cov;
        Cvar = cov * Cvar;

        % lines 47 to 49
        % generates kl biggest eigenvalues and corresponding eigenvectors
        % (captures 95% of total variance)

        fprintf('Computing eigenvalues of covariance operator ... ');
        [V,D] = eig(Cvar);
        [lambda,index] = sort(diag(D),'descend');
        V = V(:,index);
        fprintf('done.\n');
        % Use enough terms to capture 95% of fluctuation
        kl = np;
        lt = sum(lambda);
        lp = lambda(np);
        %         while lp/lt < .05,
        %             kl = kl-1;
        %             lp = lp + lambda(kl);
        %         end
        kl = problem.n;
        ratio = sum(lambda(1:kl))/sum(lambda);

        % lines 77 to 81
        % saves data for later use in truncated "KL"  expansion
        % (note the scaling of the eigenvector in line 80
        % by the square root of the associated eigenvalue ..)
%         Vl = V(:,1:kl)*diag(sqrt(lambda(1:kl)));
        Vl = V(:,1:kl)*diag((lambda(1:kl)));

        % ---------------------------------------------
        % Look at linstab_alleig_for_spinterp.m next
        %
        % lines 51 to 54
        % generate the scalar perturbation field dpsi
        % associated the multidimensional collocation point xi
        % construct the perturbed veocity field and convection matrix
        U = 1;      % normalization of velocity via max(inflow) / max(cavity b.c.)
        boundvec = ones(size(xyp(:,1)));
%         boundvec(boundp) = 0;
%         dpsi = @(xi) U*cov*(Vl*xi(:)) .* boundvec;

%         psi0 = xyp(:,1).^2 .* (1-xyp(:,2).^2) + xyp(:,2).^2;
        
        n_x = sqrt(size(xyp,1));
        xp = reshape(xyp(:,1), [n_x,n_x]);
        yp = reshape(xyp(:,2), [n_x,n_x]);
        
        Vl_xyp = reshape(Vl, [n_x,n_x,problem.n]);

%         wx_xyp = diff(Vl_xyp,2,1)./(diff(xp,2,1));
        wx_xyp = (Vl_xyp(3:end,:,:)-Vl_xyp(1:(end-2),:,:))./(xp(3:end,:)-xp(1:(end-2),:));
        % Extend on to x=-1,1 boundary
        wx_xyp_extended = zeros(size(Vl_xyp));
        wx_xyp_extended(2:end-1,:,:) = wx_xyp;
%         wx_xyp = wx_xyp(:,2:end-1,:);
%         xp = xp(2:(end-1),2:(end-1),:);
%         wy_xyp = diff(Vl_xyp,2,1)./diff(yp,2,2);
        wy_xyp = (Vl_xyp(:,3:end,:)-Vl_xyp(:,1:(end-2),:))./(yp(:,3:end)-yp(:,1:(end-2)));
        %Extend onto y=-1,1 boundary
        wy_xyp_extended = zeros(size(Vl_xyp));
        wy_xyp_extended(:,2:end-1,:) = wy_xyp;
       
%         wy_xyp = wy_xyp(2:end-1,:,:);
%         yp = yp(2:(end-1),2:(end-1),:);


        % lines 62 to 63
        % plots the perturbation field
        % Note that dpsi will be interpreted as a discrete stream function field
        % in what follows
        %
        % ---------------------------------------------
        % Look at linstab_neig_for_spinterp.m next
        %
        % lines 64 to 68
        % construct the perturbation (note multiplication by cov  in line 68)

        % line 74
        % generate the discrete convection matrix Nxi associated with the
        % perturbation (should be skew-symmetric by construction, see below)
%         psi = @(xi) psi0 + dpsi(xi);
        
        [A,N,Q,epe,eph,epw,n2] = femq1_cd_modified(xy,ev,wind_fn);
        N0 = N;

        for ii= 1:problem.n
                wind_fn = @(x,y,nel) [interp2(xp.',yp.',squeeze(wx_xyp_extended(:,:,ii)).',x,y), interp2(xp.',yp.',squeeze(wy_xyp_extended(:,:,ii)).',x,y)];
                [A,N,Q,epe,eph,epw,n2] = femq1_cd_modified(xy,ev,wind_fn);
                Nij{ii} = 0.5* ( N - N.');
        end
        Conv = @(y) N0 + combine_matrices(y, Nij); %N0 + sum(bsxfun(@times,reshape(y,1,1,[]),Nij),1);
        epsilon = 1/10;
        Diff = epsilon * A;
        AN = @(y) Diff + Conv(y); % Define the advection-diffusion matrix as a sum of fns xi of the parameters y

    otherwise
        error('Invalid test problem')
end

%% get boundary condition at current and previous times
xbd=xy(bound,1); ybd=xy(bound,2);
% bc=specific_bc(xbd,ybd);
bc=hotwallx_bc(xbd,ybd);
bc_alt=hotwallx_bc(-xbd,ybd);

invtau = 1/problem.tau;
% bc_fn = @(t,y) bc*(1-exp(-invtau*t));
bc_fn = @(t,y) bc*(1-exp(-invtau*t));
bc_fn2 =@(t,y) (t > 50).*(- bc*(1-exp(-invtau*(t-50))) + bc_alt.*(1-exp(-invtau*(t-50)))*10^(y(5)/10));

bc_fn_prime = @(t,y) invtau*bc*(exp(-invtau*t));
bc_fn2_prime = @(t,y)  (t > 50).*(- invtau*bc.*(exp(-invtau*(t-50))) + invtau*bc_alt.*(exp(-invtau*(t-50)))*10^(y(5)));

notbound = 1:size(xy,1);
notbound = setdiff(notbound, bound);

if problem.change_bc == 0
fr = @(t,y) -filter_rows_columns(AN(y),notbound,bound)*bc_fn(t,y) - Q(notbound,bound)*bc_fn_prime(t,y);
else
fr = @(t,y)  -filter_rows_columns(AN(y),notbound,bound)*bc_fn(t,y) - Q(notbound,bound)*bc_fn_prime(t,y) -filter_rows_columns(AN(y),notbound,bound)*bc_fn2(t,y) - Q(notbound,bound)*bc_fn2_prime(t,y); 
end

ANx = @(y) filter_rows_columns(AN(y),notbound,notbound);
Qx = filter_rows_columns(Q, notbound,notbound);
fx = @(t,y) fr(t,y);

fem.AN = ANx;
fem.Q = Qx;
fem.f = fx;
fem.xy = xy;
fem.ev = ev;
fem.bound = bound;
fem.notbound = notbound;
fem.bc_fn = bc_fn;

fem.mass = Q;
fem.stiffness = Diff / epsilon;
fem.conv = Conv;
[hx,hy,eex] = edgegen(xy,ev);
fem.hx = hx;
fem.hy = hy;
fem.eex = eex;
fem.ebound = ebound;
fem.lambda = epsilon;

end