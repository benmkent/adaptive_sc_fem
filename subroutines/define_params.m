function params = define_params(varargin);

params.grid = 'p1';
params.grid_param = 4; %default('Grid parameter',4);
params.grid_type = 2; %default('Grid type (2 for stretched)',2);
params.T = 100; %default('T',100);

params.timestepping = 'stabtr'; %'fixed'
params.dt0 = 1e-8;

params.t0 = -log(0.9)* 0.1; % default('t0', -log(0.9)* 0.1);
params.k_grow = 1.2;
params.k_shrink = 2;
params.k_interp = 10;
params.marking_factor = 0.1;
params.letol_lofi = 1e-1; %default('LoFi Local Error Tolerance', 1e-1);
params.dt0_lofi = params.dt0*1e1; %default('LoFi Local Error Tolerance', 1e-1);

params.letol = 1e-5; %default('Local Error Tolerance', 1e-4);
params.test_level = 8; %default('Max test level', 5);

params.rule=@(I) sum(I,2);
params.knot_fn = @(n) knots_CC(n,-1,1);
params.lev2knots = @(level) lev2knots_doubling(level);
params.l_initial = 1;
params.adapt_interp = 1;

params.downsample_to_reference = false;
params.plot = 0;
params.plottype = 'surf'; %'contourf'; % 'surf'; %
params.plotview = [-40,60];
params.plot_times = unique([logspace(-2,0,1e2), linspace(1,100,99e2+1)]);

params.jomp_estimator = 1;
params.residual_estimator = 0;

params.simplified_estimator = false;
params.reference = 0;

if nargin == 1
    params.plot = 0;
    switch varargin{1}
        case 'l4-jomp'
            params.letol = 10^(-5);
            params.dt0 = 10^(-9);
            params.adapt_interp = 1;
            params.grid = 'q1';
        case 'test-temporal-le1e0'
            params.letol = 10^(0);
            params.timestepping = 'stabtr';
        case 'test-temporal-le1e-1'
            params.letol = 10^(-1);
            params.timestepping = 'stabtr';
        case 'test-temporal-le1e-1.5'
            params.letol = 10^(-1.5);
            params.timestepping = 'stabtr';
        case 'test-temporal-le1e-2'
            params.letol = 10^(-2);
            params.timestepping = 'stabtr';
        case 'test-temporal-le1e-2.5'
            params.letol = 10^(-2.5);
            params.timestepping = 'stabtr';
        case 'test-temporal-le1e-3'
            params.letol = 10^(-3);
            params.timestepping = 'stabtr';
        case 'test-temporal-le1e-3.5'
            params.letol = 10^(-3.5);
            params.timestepping = 'stabtr';
        case 'test-temporal-le1e-4'
            params.letol = 10^(-4);
            params.timestepping = 'stabtr';
        case 'test-temporal-le1e-4.5'
            params.letol = 10^(-4.5);
            params.timestepping = 'stabtr';
        case 'test-temporal-le1e-5'
            params.letol = 10^(-5);
            params.timestepping = 'stabtr';
        case 'test-temporal-le1e-6'
            params.letol = 10^(-6);
            params.timestepping = 'stabtr';
        case 'test-temporal-ref'
            params.letol = 10^(-5);
            params.timestepping = 'stabtr';
            params.reference = 1;
        case 'test-spatial-l1'
            params.grid_param = 1;
        case 'test-spatial-l2'
            params.grid_param = 2;
        case 'test-spatial-l3'
            params.grid_param = 3;
        case 'test-spatial-l4'
            params.grid_param = 4;
        case 'test-spatial-l5'
            params.grid_param = 5;
        case 'test-spatial-l6'
            params.grid_param = 6;
        case 'test-spatial-ref'
            params.grid_param = 5;
            params.reference = 1;
        case 'test-interp-ref'
            params.reference = 1;
            params.letol = 10^(-5);
            params.grid_param = 5;
            params.l_initial = 5;
            params.k_interp = inf;
    
%             params.knot_fn = @(n) knots_leja(n,-1,1,'line');
%             params.lev2knots = @(level) lev2knots_lin(level);
        case 'test-interp'
            params.letol = 10^(-4);
            params.grid_param = 4;
            params.adapt_interp = 1;
            params.l_initial = 1;
            params.k_interp = inf;
            params.residual_estimator = 1;
%             params.knot_fn = @(n) knots_leja(n,-1,1,'line');
%             params.lev2knots = @(level) lev2knots_lin(level);
        case 'test-interp-jomp'
            params.grid_param = 6;
            params.grid_type = 1;
            params.grid = 'q1';
            params.adapt_interp = 1;
            params.l_initial = 1;
            params.residual_est = 0;
            params.marking_factor = 0.1;
            params.k_interp = 5;
            params.plot = 0;
            params.dt0=1e-2;
            params.letol = 1e-6;
            params.dt0_lofi=1e-1;
        case 'smolyak-ref'
            params.grid = 'q1';
            params.grid_param = 5; %default('Grid parameter',4);
            params.grid_type = 2; %default('Grid type (2 for stretched)',2);
            params.T = 100; %default('T',100);
            params.timestepping = 'stabtr'; %'fixed'
            params.dt0 = 1e-8;
            params.k_interp = inf;
            params.letol = 1e-7; %default('Local Error Tolerance', 1e-4);
            params.test_level = 8; %default('Max test level', 5);
            params.l_initial = 6;
            params.jomp_estimator = 1;
            params.simplified_estimator = false;
            params.reference = 1;
        case 'adaptive-base'
            params.grid = 'q1';
            params.grid_param = 5; %default('Grid parameter',4);
            params.grid_type = 2; %default('Grid type (2 for stretched)',2);
            params.T = 100; %default('T',100);
            params.timestepping = 'stabtr'; %'fixed'
            params.dt0 = 1e-8;
            params.k_interp = 1e1;
            params.letol = 1e-6; %default('Local Error Tolerance', 1e-4);
            params.test_level = 8; %default('Max test level', 5);
            params.l_initial = 1;
            params.jomp_estimator = 1;
            params.simplified_estimator = false;
            params.reference = 0;
            params.adapt_interp = 1;
        case 'smolyak-1'
            params.grid = 'q1';
            params.grid_param = 5; %default('Grid parameter',4);
            params.grid_type = 2; %default('Grid type (2 for stretched)',2);
            params.T = 100; %default('T',100);
            params.timestepping = 'stabtr'; %'fixed'
            params.dt0 = 1e-8;
            params.k_interp = inf;
            params.letol = 1e-6; %default('Local Error Tolerance', 1e-4);
            params.test_level = 8; %default('Max test level', 5);
            params.l_initial = 1;
            params.jomp_estimator = 1;
            params.simplified_estimator = false;
            params.reference = 0;
        case 'smolyak-2'
            params.grid = 'q1';
            params.grid_param = 5; %default('Grid parameter',4);
            params.grid_type = 2; %default('Grid type (2 for stretched)',2);
            params.T = 100; %default('T',100);
            params.timestepping = 'stabtr'; %'fixed'
            params.dt0 = 1e-8;
            params.k_interp = inf;
            params.letol = 1e-6; %default('Local Error Tolerance', 1e-4);
            params.test_level = 8; %default('Max test level', 5);
            params.l_initial = 2;
            params.jomp_estimator = 1;
            params.simplified_estimator = false;
            params.reference = 0;
        case 'smolyak-3'
            params.grid = 'q1';
            params.grid_param = 5; %default('Grid parameter',4);
            params.grid_type = 2; %default('Grid type (2 for stretched)',2);
            params.T = 100; %default('T',100);
            params.timestepping = 'stabtr'; %'fixed'
            params.dt0 = 1e-8;
            params.k_interp = inf;
            params.letol = 1e-6; %default('Local Error Tolerance', 1e-4);
            params.test_level = 8; %default('Max test level', 5);
            params.l_initial = 3;
            params.jomp_estimator = 1;
            params.simplified_estimator = false;
            params.reference = 0;
        case 'smolyak-4'
            params.grid = 'q1';
            params.grid_param = 5; %default('Grid parameter',4);
            params.grid_type = 2; %default('Grid type (2 for stretched)',2);
            params.T = 100; %default('T',100);
            params.timestepping = 'stabtr'; %'fixed'
            params.dt0 = 1e-8;
            params.k_interp = inf;
            params.letol = 1e-6; %default('Local Error Tolerance', 1e-4);
            params.test_level = 8; %default('Max test level', 5);
            params.l_initial = 4;
            params.jomp_estimator = 1;
            params.simplified_estimator = false;
            params.reference = 0;
        case 'smolyak-5'
            params.grid = 'q1';
            params.grid_param = 5; %default('Grid parameter',4);
            params.grid_type = 2; %default('Grid type (2 for stretched)',2);
            params.T = 100; %default('T',100);
            params.timestepping = 'stabtr'; %'fixed'
            params.dt0 = 1e-8;
            params.k_interp = inf;
            params.letol = 1e-6; %default('Local Error Tolerance', 1e-4);
            params.test_level = 8; %default('Max test level', 5);
            params.l_initial = 5;
            params.jomp_estimator = 1;
            params.simplified_estimator = false;
            params.reference = 0;
    end
end
