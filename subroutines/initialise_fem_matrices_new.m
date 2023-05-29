function fem = initialise_fem_matrices_new(problem, params)
%% Define geometry of spatial domain

%% Generate grid and mass matrix
switch params.grid
    case 'q1'
        pde=2; domain=1;
        [mv, xy, bound, mbound, grid_type, outbc, x, y] = ...
            rsquare_domain_modified(params.grid_param, params.grid_type);

        [ev,ebound] = q1grid(xy,mv,bound,mbound);
        [Q] = femq1_mass(xy,ev);
        [hx,hy,eex] = edgegen(xy,ev);
        tve=nan;
    case 'p1'
        %         square_domain(2,1) ;
        % initial grid param
        ini_grid_param = 2;
        [ mv, xy, bound, mbound, grid_type, outbc, x, y] = t_square_domain_modified(2, ini_grid_param, params.grid_type);
        [ev,ebound] = p1grid(xy,mv,bound,mbound);
        for i=1:(params.grid_param-1),
            [xy,ev,bound,ebound] = p1_refinement(xy,ev,bound,ebound);
        end
        % validate mesh
        [Q] = femp1_mass(xy,ev);
        [eex,tve,hx] = tedgegen(xy,ev);
        hy = 0*hx;

        T = ev;
        %         T = delaunay(xy(:,1),xy(:,2));
        %         T = delaunayTriangulation(xy(:,1),xy(:,2));
        fem.T = T;
    otherwise
        error('Only P1 or Q1');
end

%% Preprocess
if ~problem.linear
    if strcmp(problem.wind_fn, 'groundwater') & strcmp(problem.diff_fn, 'groundwater')
        %% Set up solver for Darcy problem to give parameter dependent solution for flow field
        problem.u = @(x,y) darcy_flow_sln(x,y, fem);

    end
end

%% Define wind field and matrices
if problem.linear
    switch problem.wind_fn
        case 'zero'
            wind_fn{1} = @(x,y,nel) ones(size(x,1),2) * 0;
            kk=1;
        case 'recirculating'
            wind_fn{1} = @(x,y,nel) wind_fn_parameterised(x,y,[0,0], 1);
            kk=1;
        case 'recirculating-alt'
            wind_fn{1} = @(x,y,nel) wind_fn_parameterised(x,y,[0.5,0.5], 2);
            kk=1;
        case 'recirculating-uncertain'
            n_dim = problem.nWind;
            sigma = problem.sigmaWind;
            sqrt_n = round(sqrt(n_dim));
            dist_centres = 2/sqrt_n;
            centres = -1 + dist_centres*0.5 + dist_centres * (0:(sqrt_n-1));

            wind_fn{1} = @(x,y,nel) wind_fn_parameterised(x,y,[0,0], 1);

            kk=1;
            for ii= 1:sqrt_n
                for jj = 1:sqrt_n
                    kk=kk+1;
                    x_w = centres(ii);
                    y_w = centres(jj);
                    wind_fn{kk} = @(x,y,nel) sigma* wind_fn_parameterised(x,y,[x_w,y_w], sqrt_n);
                end
            end
        case 'recirculating-uncertain-decay'
            n_dim = problem.nWind;
            sigma = problem.sigmaWind;
            sqrt_n = round(sqrt(n_dim));
            dist_centres = 2/sqrt_n;
            centres = -1 + dist_centres*0.5 + dist_centres * (0:(sqrt_n-1));

            wind_fn{1} = @(x,y,nel) wind_fn_parameterised(x,y,[0,0], 1);

            kk=1;
            for ii= 1:sqrt_n
                for jj = 1:sqrt_n
                    kk=kk+1;
                    x_w = centres(ii);
                    y_w = centres(jj);
                    wind_fn{kk} = @(x,y,nel) sigma*0.7^((n_dim-kk+2)) * wind_fn_parameterised(x,y,[x_w,y_w], sqrt_n);
                end
            end
        case 'jomp'
            n_dim = problem.nWind;
            sigma = problem.sigmaWind;

            kk=1;
            wind_fn{1} = @(x,y,nel) wind_fn_parameterised(x,y,[0,0], 1);

            % Set up four large recirculating, then 16 small
            sqrt_n = 2; %round(sqrt(n_dim));
            dist_centres = 2/sqrt_n;
            centres = -1 + dist_centres*0.5 + dist_centres * (0:(sqrt_n-1));

            for ii= 1:sqrt_n
                for jj = 1:sqrt_n
                    kk=kk+1;
                    x_w = centres(ii);
                    y_w = centres(jj);
                    wind_fn{kk} = @(x,y,nel) sigma * wind_fn_parameterised(x,y,[x_w,y_w], sqrt_n);
                end
            end
            % now 16
            sqrt_n = 4; %round(sqrt(n_dim));
            dist_centres = 2/sqrt_n;
            centres = -1 + dist_centres*0.5 + dist_centres * (0:(sqrt_n-1));

            % Create anisotropic coeffs and permute to (seeded) random order
            coeffs_ani =  (1.2).^(-1:-1:-16);
            s = rng(0);
            r = randperm(16);
            coeffs_ani = coeffs_ani(r);
            kk2 = 0;
            for ii= 1:sqrt_n
                for jj = 1:sqrt_n
                    kk=kk+1;
                    kk2 = kk2+1;
                    x_w = centres(ii);
                    y_w = centres(jj);
                    wind_fn{kk} = @(x,y,nel) sigma * coeffs_ani(kk2) * wind_fn_parameterised(x,y,[x_w,y_w], sqrt_n);
                end
            end
        case 'jomp32'
            n_dim = problem.nWind;
            sigma = problem.sigmaWind;
            %             sqrt_n = round(sqrt(n_dim));
            dist_centres = 2/8;
            centres_x = -1 + dist_centres*0.5 + dist_centres * (0:(8-1));
            centres_x(2:2:end) = [];
            dist_centres = 2/8;
            centres_y = -1 + dist_centres*0.5 + dist_centres * (0:(8-1));

            wind_fn{1} = @(x,y,nel) wind_fn_parameterised(x,y,[0,0], 1);

            coeffs([3	28	20	7	2	8	23	4])  = 1e0;
            coeffs([16	24	6	18	19	14	29	17]) = 1e-1;
            coeffs([22	26	12	1	25	31	30	32]) = 1e-2;
            coeffs([13	5	9	21	15	11	10	27]) = 1e-3;

            kk=1;
            for ii= 1:4
                for jj = 1:8
                    kk=kk+1;
                    x_w = centres_x(ii);
                    y_w = centres_y(jj);
                    wind_fn{kk} = @(x,y,nel) coeffs(kk-1)* wind_fn_parameterised(x,y,[x_w,y_w], 8);
                end
            end
        case 'jomp64'
            n_dim = problem.nWind;
            sigma = problem.sigmaWind;
            sqrt_n = round(sqrt(n_dim));
            dist_centres = 2/sqrt_n;
            centres = -1 + dist_centres*0.5 + dist_centres * (0:(sqrt_n-1));

            wind_fn{1} = @(x,y,nel) wind_fn_parameterised(x,y,[0,0], 1);

            coeffs([51	44	43	64	7	48	20	21	53	6	57	29	12	26	30	24])  = 1e0;
            coeffs([14	22	59	3	38	55	52	10	28	33	41	63	50	32	8	15]) = 1e-1;
            coeffs([56	45	19	36	42	49	35	11	58	31	60	5	4	61	16	1]) = 1e-2;
            coeffs([62	13	40	2	39	46	54	25	34	23	17	37	27	47	18	9]) = 1e-3;
            kk=1;
            for ii= 1:sqrt_n
                for jj = 1:sqrt_n
                    kk=kk+1;
                    x_w = centres(ii);
                    y_w = centres(jj);
                    wind_fn{kk} = @(x,y,nel) coeffs(kk-1)* wind_fn_parameterised(x,y,[x_w,y_w], sqrt_n);
                end
            end
        case 'djs-rf'
            n_dim = problem.nWind;
            cov = 5e0;

            wind_fn{1} = @(x,y,nel) wind_fn_parameterised(x,y,[0,0], 1);

            [wind_fn(2:(n_dim+1)), lambda, Vl, w_max, w_mean, xp, yp] = djs_wind_perturbation(n_dim,cov,7);

            save('rf-construction','wind_fn','lambda','Vl','w_max','w_mean','xp','yp');

            kk = n_dim+1;
        case 'xonly'
            wind_fn{1} = @(x,y,nel) [1,0];
            kk=1;
        case 'xonly-weak'
            xLength = 1;
            xTime = 1;
            xU = xLength/xTime;
            wind_fn{1} = @(x,y,nel) [xU,0];
            kk=1;
        case 'xonly-perturbed'
            xLength = 1;
            xTime = 20;
            theta = 0;
            U = xLength/xTime;
            xU = U *cos(theta);
            yU = U *sin(theta);
            wind_fn{1} = @(x,y,nel) [xU,yU];
            kk=1;

            n_dim = problem.nWind;
            sigma = problem.sigmaWind;
            sqrt_n = round(sqrt(n_dim));
            dist_centres = 2/sqrt_n;
            centres = -1 + dist_centres*0.5 + dist_centres * (0:(sqrt_n-1));
            for ii= 1:sqrt_n
                for jj = 1:sqrt_n
                    kk=kk+1;
                    x_w = centres(ii);
                    y_w = centres(jj);
                    wind_fn{kk} = @(x,y,nel) sigma  * U * wind_fn_parameterised(x,y,[x_w,y_w], sqrt_n);
                end
            end
        case 'pouiseuille'
            wind_fn{1} = @(x,y,nel) 0.1*[(1-y.^2) , 0*y];
            kk=1;
        otherwise
            error('Unknown wind fn');
    end
    % Append extra zero wind functions as necessary to get correct number of rv
    if kk <= problem.nWind
        for kk_ii = (kk+1):problem.nWind+1
            wind_fn{kk_ii} = @(x,y,nel) ones(size(x,1),2) * 0;
        end
    end

    % plot_wind_fn(wind_fn);

    for ii = 1:length(wind_fn)
        fprintf('(%d of %d)...',ii,length(wind_fn));
        switch params.grid
            case 'q1'
                N{ii} = femq1_conv(xy,ev, wind_fn{ii});
            case 'p1'
                N{ii} = femp1_conv(xy,ev, wind_fn{ii});
        end
    end

    nYConv = length(wind_fn)-1;

    if nYConv > 0
        Conv = @(y) N{1} + combine_matrices(y(1:nYConv), N(2:(nYConv+1)));
    else
        Conv = @(y) N{1};
    end

    %     wind_fn_eval = @(x1,x2, y)  sum ( [1; y(1:nYConv)] .* cell2mat(cellfun(@(fi)fi(x1,x2),wind_fn,'UniformOutput',false)'),1);
    wind_fn_eval = @(x1,x2, y_in) evaluate_w_lincomb(x1,x2,[1;y_in((1):(nYConv))],wind_fn);

else
    switch problem.wind_fn
        case 'groundwater'
            u = problem.u;

            porosity = 1;

            wind_fn = @(x1,x2,y) u([x1,x2],z) / porosity;
    end

    switch params.grid
        case 'q1'
            Conv = @(y) femq1_conv(xy,ev, @(x1,x2) wind_fn(x1,x2,y));
        case 'p1'
            Conv = @(y) femp1_conv(xy,ev, @(x1,x2) wind_fn(x1,x2,y));
    end

    wind_fn_eval = @(x1,x2, y)  wind_fn(x1,x2,y);
end

%% Set up diffusion
if problem.linear == 1
    switch problem.diff_fn
        case 'fixed'
            a_fn{1} = @(x,y) problem.viscosity * [ones(size(x,1),1), zeros(size(x,1),2), ones(size(x,1),1)];
            kk = 1;
        case 'scalar-var'
            a_fn{1} = @(x,y) problem.viscosity * [ones(size(x,1),1), zeros(size(x,1),2), ones(size(x,1),1)];
            a_fn{2} = @(x,y) problem.viscosity_var * [ones(size(x,1),1), zeros(size(x,1),2), ones(size(x,1),1)];
            kk = 2;
        case 'eigel'
            nDim = problem.nDiff;

            alpha_bar = 0.547;
            k = @(mm) floor(-0.5 + sqrt(1/4 + 2*mm));
            beta_1 = @(mm) mm - k(mm)*(k(mm)+1)/2;
            beta_2 = @(mm) k(mm) - beta_1(mm);
            alpha_m = @(mm) alpha_bar * mm^-2;

            a_fn{1} = @(x,y) eye(2);
            xmap = @(x) (x+1)/2;
            for mm = 1:nDim
                diff_fn = @(x1,x2) alpha_m(mm) .* cos(2*pi* beta_1(mm) *x1) .* cos(2*pi*beta_2(mm)*x2);
                diff_fn_shifted = @(x1,x2) diff_fn(xmap(x1),xmap(x2))  * eye(2);
                a_fn{1+mm} = diff_fn_shifted;
            end
            kk=1+mm;

        case 'eigel-y'
            nDim = problem.nDiff;

            alpha_bar = 0.547;
            k = @(mm) floor(-0.5 + sqrt(1/4 + 2*mm));
            beta_1 = @(mm) mm - k(mm)*(k(mm)+1)/2;
            beta_2 = @(mm) k(mm) - beta_1(mm);
            alpha_m = @(mm) alpha_bar * mm^-2;

            a_fn{1} = @(x,y) eye(2);
            xmap = @(x) (x+1)/2;
            for mm = 1:nDim
                diff_fn = @(x1,x2) alpha_m(mm) .* cos(2*pi* beta_1(mm) *x1) .* cos(2*pi*beta_2(mm)*x2);
                diff_fn_shifted = @(x1,x2) [alpha_m(mm), 0; 0, diff_fn(xmap(x1),xmap(x2)) ];
                a_fn{1+mm} = diff_fn_shifted;
            end
            kk=1+mm;

        case 'vidlickova'
            nYDiff = problem.nDiff;
            %             if nYDiff ~= 2; error('Problem set up wrong'); end

            a0 = 1;
            a_fn{1} = @(x,y) a0* eye(2);
            xmap = @(x) (x+1)/2;

            for mm = 1:2
                a_fn{1+mm} = @(x,y) eye(2) * (cos(2*pi*mm*xmap(x)) + cos(2*pi*mm*xmap(y)))/(pi*mm)^2;
            end
            kk=3;

        case 'complex'
            a0 = 1e-4;
            a_fn{1} = @(x,y) a0 * [1,0, 0, 1i];
    end

    % Append extra zero diff functions as necessary to get correct number of rv
    if kk <= problem.nDiff
        for kk_ii = (kk+1):problem.nDiff+1
            a_fn{kk_ii} = @(x,y,nel) 0 * [ones(size(x,1),1), zeros(size(x,1),2), ones(size(x,1),1)];
        end
    end

    for ii = 1:length(a_fn)
        fprintf('(%d of %d)...',ii,length(a_fn));
        switch params.grid
            case 'q1'
                [A{ii}, H] = femq1_diff_only(xy,ev, a_fn{ii});
            case 'p1'
                [A{ii}, H] = femp1_diff_only(xy,ev, a_fn{ii});
        end
    end

    nYDiff = length(a_fn)-1;

    if nYDiff > 0
        Diff = @(y) A{1} + combine_matrices(y((nYConv+1):(nYConv+nYDiff)), A(2:end));
    else
        Diff = @(y) A{1};
    end

    %     a_fn_eval = @(x1,x2, y)  [1,0] * reshape(sum ( [1; y((nYConv+1):(nYConv+nYDiff))] .* cell2mat(cellfun(@(fi)reshape(fi(x1,x2),[1,4]),a_fn,'UniformOutput',false)'),1),[2,2]) *[1;0];
    a_fn_eval = @(x1,x2, y) evaluate_a_lincomb(x1,x2,[1;y((nYConv+1):(nYConv+nYDiff))],a_fn);

    % Construct pairwise integration matrices
    %     switch params.grid
    %         case 'q1'
    %             [a_pairs] = femq1_diff_pairs(xy,ev, a_fn);
    %         case 'p1'
    %             [a_pairs] = femp1_diff_pairs(xy,ev, a_fn);
    %     end
    % a_pairs = [];

else
    switch problem.diff_fn
        case 'groundwater'
            u = problem.u;

            a_fn = @(x1,x2,y) (aL-aT) * u([x1,x2],y)' * u([x1,x2],z) / norm(u([x1,x2],y),2) +...
                (aT*norm(u([x1,x2],y),2) + Dd) * ones(2,2);
    end
    switch params.grid
        case 'q1'
            [~, H] = femq1_diff_only(xy,ev, @(x1,x2) a_fn(x1,x2,zeros(problem.n,1)));
            Diff = @(y) femq1_diff_only(xy,ev, @(x1,x2) a_fn(x1,x2,y));
        case 'p1'
            [~, H] = femp1_diff_only(xy,ev, @(x1,x2) a_fn(x1,x2,zeros(problem.n,1)));
            Diff = @(y) femp1_diff_only(xy,ev, @(x1,x2) a_fn(x1,x2,y));
    end

    a_fn_eval = @(x1,x2,y) a_fn(x1,x2,y);
end


% n_test_lambda = 1e2;
% x_lambda = cos(linspace(0,2*pi,n_test_lambda));
% y_lambda = sin(linspace(0,2*pi,n_test_lambda));
% for ii = 1:n_test_lambda
%     v_lambda = [x_lambda(ii);y_lambda(ii)];
%     lambda_test(ii) = [v_lambda' * a()]

% lambda = problem.viscosity;
lambda = 0.1;

%% Set up boundary condition
notbound = 1:size(xy,1);
notbound = setdiff(notbound, bound);

xbd=xy(bound,1); ybd=xy(bound,2);
switch problem.bc
    case 'dirichlet'
        bc= zeros(size(xbd,1),1);
        bc_fn = @(t,y) bc*zeros(1,length(t));
        bc_fn_prime = @(t,y) bc*zeros(1,length(t));
    case 'hotwall'
        bc=hotwallx_bc(xbd,ybd);
        invtau = 1/problem.tau;
        % bc_fn = @(t,y) bc*(1-exp(-invtau*t));
        bc_fn = @(t,y) bc*(1-exp(-invtau*t));
        bc_fn_prime = @(t,y) invtau*bc*(exp(-invtau*t));
    case 'hotwall-timebc'
        bc=hotwallx_bc(xbd,ybd);
        invtau = 1/problem.tau;
        % bc_fn = @(t,y) bc*(1-exp(-invtau*t));
        bc_fn = @(t,y) bc*((1-exp(-invtau*t)).*(1+0.2*sin(2*pi*t)));
        bc_fn_prime = @(t,y) invtau*bc*((exp(-invtau*t)).*(1+0.2*sin(2*pi*t))) + 2*pi*bc*((1-exp(-invtau*t)).*(0.2*cos(2*pi*t)));        
    case 'neumann'
        bound = find(abs((xy(:,2) - -1)) < 1e-8);
        notbound = 1:size(xy,1);
        notbound = setdiff(notbound, bound);
        xbd=xy(bound,1); ybd=xy(bound,2);
        bc_fn = @(t,y) zeros(length(bound),length(t));
        bc_fn_prime =  @(t,y) zeros(length(bound),length(t));
    case 'pipe'
        bound = sort([find(abs((xy(:,2) - -1)) < 1e-8); find(abs((xy(:,2) - 1)) < 1e-8)]);
        notbound = 1:size(xy,1);
        notbound = setdiff(notbound, bound);
        xbd=xy(bound,1); ybd=xy(bound,2);
        bc_fn = @(t,y) zeros(length(bound),length(t));
        bc_fn_prime =  @(t,y) zeros(length(bound),length(t));
    otherwise
        error('BC does not exist');
end

%% Set up forcing function
switch problem.f
    case 'zero'
        force_fn = @(x) 0; %zeros(size(xy,1),1);
    case 'square'
        x1_source = -0.5;
        x2_source = 0;
        r = 0.05;
        f_factor = 20;

        force_fn = @(x) f_factor * double(abs(x(:,1)-x1_source) < r & abs(x(:,2)-x2_source) < r);
    case 'circle'
        x1_source = -0.5;
        x2_source = 0;
        r = 0.05;
        f_factor = 10;

        force_fn = @(x) f_factor*double((x(:,1)-x1_source).^2 + (x(:,2)-x2_source).^2 < r^2);
    case 'one'
        force_fn = @(x) ones(size(x,1),1);
    case 'vidlickova'
        xmap = @(x) (x+1)/2;

        x1_source = 0.5; % [0,1]x[0,1] coords
        x2_source = 0.5;
        r = 0.1;

        f_factor = 20;
        force_fn = @(x) f_factor * double(abs(xmap(x(:,1))-x1_source) < r & abs(xmap(x(:,2))-x2_source) < r);
    otherwise
        error('Unknown forcing')
end

switch params.grid
    case 'p1'
        f = femp1_force(xy,ev, force_fn);
    case 'q1'
        f = femq1_force(xy,ev, force_fn);
end

%% Create fem structure
fem.wind_fn = wind_fn_eval;
fem.a_fn = a_fn_eval;
% fem.a_pairs = a_pairs;

fem.diff = Diff;
fem.conv = Conv;
fem.mass = Q;
fem.stiffness = H;
fem.bc_fn = bc_fn;
fem.bc_fn_prime = bc_fn_prime;
fem.f = f;
fem.lambda = lambda;

fem.xy = xy;
fem.ev = ev;
fem.bound = bound;
fem.notbound = notbound;
fem.hx = hx;
fem.hy = hy;
fem.eex = eex;
fem.ebound = ebound;
fem.tve = tve;

fprintf('Spatial DOF per pt %d\n', size(fem.xy,1));

end
