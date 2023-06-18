%% Define problem
problem = define_problem('doubleglazing');

%% Run experiments

reference = define_reference('none');
params = define_params('l4-jomp');

params.adapt_interp = 0; % Interpolate for new collocation points
params.k_interp = 10; % c_{tol} safety factor
params.marking_factor = 0.1; % Marking factor
params.letol = 1e-4; % Local error tolerance

adaptive_sc_fem;
save(['l4-jomp.mat'],'reference','data_table','fem','problem','params','-v7.3')

plot_data(data_table, reference,fem,params,problem)