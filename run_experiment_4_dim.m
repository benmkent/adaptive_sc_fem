%% Define problem
problem = define_problem('doubleglazing');

%% Set up reference
params = define_params('l4-jomp');
params.l_initial = 5;
params.letol = 1e-7;
adaptive_sc_fem;
save(['reference.mat'],'reference','data_table','fem','problem','params', '-v7.3')

%% Run experiments

reference = define_reference('test-folder');
params = define_params('l4-jomp');

params.adapt_interp = 1; % Interpolate for new collocation points
params.simplified_estimator = 1; % Alternative error estimator (all colloc pts use the same)
params.k_interp = 10; % c_{tol} safety factor
params.marking_factor = 0.1; % Marking factor
params.letol = 1e-4; % Local error tolerance

adaptive_sc_fem;
save(['l4-jomp.mat'],'reference','data_table','fem','problem','params','-v7.3')

plot_data(data_table, reference,fem,params,problem)