function fem = initialise_fem_matrices(problem, params,fem)
%INITIALISE_FEM_MATRICES For a FEM mesh construct discrete adv-diff problem
% Inputs
%   problem     defines the adv-diff problem
%   params      defines the approximation parameters
%   fem         contains the FEM mesh properties
%
% Output
%   fem         updated fem structure

%% Set up FEM mesh variables
xy = fem.xy;
ev = fem.ev;
bound = fem.bound;
notbound = fem.notbound;

%% Preprocess

%% Define wind field and matrices
if problem.linear
    %% Define fields
    switch problem.wind_fn
        case 'zero'
            wind_fn{1} = @(x,y,nel) ones(size(x,1),2) * 0;
            kk=1;
        case 'recirculating'
            wind_fn{1} = @(x,y,nel) recirculating_parameterised(x,y,[0,0], 1);
            kk=1;
        case 'recirculating-alt'
            wind_fn{1} = @(x,y,nel) recirculating_parameterised(x,y,[0.5,0.5], 2);
            kk=1;
        case 'recirculating-uncertain'
            n_dim = problem.nWind;
            sigma = problem.sigmaWind;
            % Set up parameters for offset "eddie" type wind field
            sqrt_n = round(sqrt(n_dim));
            dist_centres = 2/sqrt_n;
            centres = -1 + dist_centres*0.5 + dist_centres * (0:(sqrt_n-1));

            % Define mean field
            wind_fn{1} = @(x,y,nel) recirculating_parameterised(x,y,[0,0], 1);

            % For each parameter define a smaller circulating wind field
            % with magnitude scaled by $\sigma$
            kk=1;
            for ii= 1:sqrt_n
                for jj = 1:sqrt_n
                    kk=kk+1;
                    x_w = centres(ii);
                    y_w = centres(jj);
                    wind_fn{kk} = @(x,y,nel) sigma* recirculating_parameterised(x,y,[x_w,y_w], sqrt_n);
                end
            end
        case 'djs-rf'
            % Construct an advection field using an approximated stream
            % function
            n_dim = problem.nWind;
            cov = 1e0;

            % Define mean recirculating wind field.
            wind_fn{1} = @(x,y,nel) wind_fn_parameterised(x,y,[0,0], 1);

            % Define perturbing fields through stream function
            [wind_fn(2:(n_dim+1)), lambda, Vl, w_max, w_mean, xp, yp] = djs_wind_perturbation(n_dim,cov,7);
            %             save('rf-construction','wind_fn','lambda','Vl','w_max','w_mean','xp','yp');

            kk = n_dim+1;
        case 'xonly'
            % Construct an advection field that is w(x,y)=[1,0]
            wind_fn{1} = @(x,y,nel) [1,0];
            kk=1;
        case 'xonly-weak'
            % Construct an advection field that is w(x,y)=[0.1,0]
            xLength = 1;
            xTime = 10;
            xU = xLength/xTime;
            wind_fn{1} = @(x,y,nel) [xU,0];
            kk=1;
        case 'xonly-perturbed'
            % Construct an advection field that is rotated w(x,y)=[0.1,0]
            xLength = 1;
            xTime = 20;
            theta = 0;
            U = xLength/xTime;
            xU = U *cos(theta);
            yU = U *sin(theta);
            wind_fn{1} = @(x,y,nel) [xU,yU];
            kk=1;

            % perturb wind field with recirculating subfields
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
                    wind_fn{kk} = @(x,y,nel) sigma  * U * recirculating_parameterised(x,y,[x_w,y_w], sqrt_n);
                end
            end
        case 'pouiseuille'
            % Construct pouiseille flow no parameters.
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

    %% Construct the advection matrices
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

    %% Construct affine field
    if nYConv > 0
        Conv = @(y) N{1} + combine_matrices(y(1:nYConv), N(2:(nYConv+1)));
    else
        Conv = @(y) N{1};
    end

    wind_fn_eval = @(x1,x2, y_in) evaluate_w_lincomb(x1,x2,[1;y_in((1):(nYConv))],wind_fn);

else
    %% Non affine advection
    % Can define through a predefined non-affine function u
    % BMK - Should change this notation
    u = problem.u;
    wind_fn = @(x1,x2,y) u([x1,x2],z);

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
    % Sets up diffusion field as a Cartesian tensor.
    % Note -- this is not fully implemented! For error estimation in
    % particular the field should be assumed to be isotropic, with zeros on
    % the off-diagonals (index 2 and 3) and equal values in indices 1 and 2.
    switch problem.diff_fn
        case 'fixed'
            % Fixed diffusion coefficient.
            a_fn{1} = @(x,y) problem.viscosity * [ones(size(x,1),1), zeros(size(x,1),2), ones(size(x,1),1)];
            kk = 1;
        case 'scalar-var'
            % Scalar diffusion with one parameter.
            a_fn{1} = @(x,y) problem.viscosity * [ones(size(x,1),1), zeros(size(x,1),2), ones(size(x,1),1)];
            a_fn{2} = @(x,y) problem.viscosity_var * [ones(size(x,1),1), zeros(size(x,1),2), ones(size(x,1),1)];
            kk = 2;
        case 'eigel'
            % WARNING - Not been tested properly.
            % This sets up the so called Eigel diffusion field from
            % https://doi.org/10.1016/j.cma.2013.11.015
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

        case 'complex'
            % Warning - this has not been tested properly.
            a0 = 1e-4;
            a_fn{1} = @(x,y) a0 * [1,0, 0, 1i];
    end

    % Append extra zero diff functions as necessary to get correct number of rv
    if kk <= problem.nDiff
        for kk_ii = (kk+1):problem.nDiff+1
            a_fn{kk_ii} = @(x,y,nel) 0 * [ones(size(x,1),1), zeros(size(x,1),2), ones(size(x,1),1)];
        end
    end

    %% Compute diffusion matrices
    for ii = 1:length(a_fn)
        fprintf('(%d of %d)...',ii,length(a_fn));
        switch params.grid
            case 'q1'
                [A{ii}, H] = femq1_diff_only(xy,ev, a_fn{ii});
            case 'p1'
                [A{ii}, H] = femp1_diff_only(xy,ev, a_fn{ii});
        end
    end

    %% Construct the affine diffusion matrix.
    nYDiff = length(a_fn)-1;

    if nYDiff > 0
        Diff = @(y) A{1} + combine_matrices(y((nYConv+1):(nYConv+nYDiff)), A(2:end));
    else
        Diff = @(y) A{1};
    end

    %     a_fn_eval = @(x1,x2, y)  [1,0] * reshape(sum ( [1; y((nYConv+1):(nYConv+nYDiff))] .* cell2mat(cellfun(@(fi)reshape(fi(x1,x2),[1,4]),a_fn,'UniformOutput',false)'),1),[2,2]) *[1;0];
    a_fn_eval = @(x1,x2, y) evaluate_a_lincomb(x1,x2,[1;y((nYConv+1):(nYConv+nYDiff))],a_fn);

else
    %% Construct diffusion field for each y by explicitly calling the appropriate FEM function.
    u = problem.u;
    a_fn = @(x1,x2,y) (aL-aT) * u([x1,x2],y)' * u([x1,x2],z) / norm(u([x1,x2],y),2) +...
        (aT*norm(u([x1,x2],y),2) + Dd) * ones(2,2);
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

%% Set up boundary condition
% Identify boundary indices and coordinates
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
        % Define homgenous Dirichlet BC at x_2=-1.
        bound = find(abs((xy(:,2) - -1)) < 1e-8);
        notbound = 1:size(xy,1);
        notbound = setdiff(notbound, bound);
        xbd=xy(bound,1); ybd=xy(bound,2);
        bc_fn = @(t,y) zeros(length(bound),length(t));
        bc_fn_prime =  @(t,y) zeros(length(bound),length(t));
    case 'pipe'
        % Define homegenous Dirichlet BC at x_2=+-1
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
    otherwise
        error('Unknown forcing')
end

%% Construct forcing vector 
switch params.grid
    case 'p1'
        f = femp1_force(xy,ev, force_fn);
    case 'q1'
        f = femq1_force(xy,ev, force_fn);
end

%% Update fem structure
% Store constructed matrices and functions into FEM structure.
fem.wind_fn = wind_fn_eval;
fem.a_fn = a_fn_eval;

fem.diff = Diff;
fem.conv = Conv;
fem.mass = fem.Q;
fem.stiffness = H;
fem.bc_fn = bc_fn;
fem.bc_fn_prime = bc_fn_prime;
fem.f = f;

fprintf('Spatial DOF per pt %d\n', size(fem.xy,1));

end
