clear all;

addpath(genpath('sparse-grids-matlab-kit_v-18-10_Esperanza'));

% problem.name ='one-d'; problem.n = 1;
% problem.name ='quadrants-isotropic'; problem.n=4;
% problem.name ='quadrants-anisotropic'; problem.n=4;
% problem.name='djs-div-free';  problem.n=3; problem.corrx=0.1; problem.corry=1; problem.delta=1; problem.cov=0.5;
problem.name = 'eddies'; problem.n = 4; problem.sigma=0.5;

params.lte_tol = 1e-4;
params.solvertype = 'stabtr';
params.dt0 = 1e-9;
params.test_level = 10;

params.tau = 1e-3;
params.hot_pc = 0.9;
params.T = 100;
params.nTimes = 40;
params.grid_param=4;
params.grid_type=2;

params.compute_interp_norms = true;

params.ge_tol = 1e-6;
params.interp_tol = 1e-5;
params.rule=@(I) sum(I,2);
params.l_initial = 4;
params.adapt_type = 'zero';
C=multiidx_gen(problem.n,params.rule,problem.n + params.l_initial-1,1);
params.C_initial = C;

%% Run First Pass For Reference Approximation
results_file=[mfilename '_reference'];
run_sc_functions;

%% Run for varied C
params.ge_tol = 1e-5;

results_file=[mfilename  '_AdaptZero_1em4'];
params.l_initial = 1;
params.interp_tol = 1e-4;
params.adapt_type = 'zero';
C=multiidx_gen(problem.n,params.rule,problem.n + params.l_initial-1,1);
params.C_initial = C;
run_sc_functions;

results_file=[mfilename  '_AdaptInterp_1em4'];
params.l_initial = 1;
params.interp_tol = 1e-4;
params.adapt_type = 'interp-all';
C=multiidx_gen(problem.n,params.rule,problem.n + params.l_initial-1,1);
params.C_initial = C;
run_sc_functions;

params.interp_tol = inf;

results_file=[mfilename  '_Linit1'];
params.l_initial = 1;
C=multiidx_gen(problem.n,params.rule,problem.n + params.l_initial-1,1);
params.C_initial = C;
run_sc_functions;

results_file=[mfilename  '_Linit2'];
params.l_initial = 2;
C=multiidx_gen(problem.n,params.rule,problem.n + params.l_initial-1,1);
params.C_initial = C;
run_sc_functions;

results_file=[mfilename  '_Linit3'];
params.l_initial = 3;
C=multiidx_gen(problem.n,params.rule,problem.n + params.l_initial-1,1);
params.C_initial = C;
run_sc_functions;

results_file=[mfilename  '_Linit4'];
params.l_initial = 4;
C=multiidx_gen(problem.n,params.rule,problem.n + params.l_initial-1,1);
params.C_initial = C;
run_sc_functions;

zip(['results/zipped/' mfilename '.zip'], {['results/' mfilename '*.mat']})
delete(['results/' mfilename '*.mat']);
