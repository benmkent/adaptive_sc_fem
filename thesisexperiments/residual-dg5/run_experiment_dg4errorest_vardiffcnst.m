addpath ../../adaptive_sc_fem/;
problem = define_problem('doubleglazing-vardiff');

%% Set up reference
params = define_params('smolyak-ref');

params.grid = 'p1';
params.grid_param = 5;
params.grid_type=1;
params.jomp_estimator = 1;
params.residual_estimator = 1;

adaptive_sc_fem;
save(['reference.mat'],'reference','data_table','fem','problem','params', '-v7.3')

%% Run experiments
names = {'smolyak-1','smolyak-2','smolyak-3','smolyak-4','smolyak-5'};

for name = names
    reference = define_reference('test-folder');
    params = define_params(name{1});

    params.grid = 'p1';
    params.grid_param = 5;
    params.grid_type=1;
    params.jomp_estimator = 1;
    params.residual_estimator = 1;


    adaptive_sc_fem;
    save([name{1} '.mat'],'reference','data_table','fem','problem','params','-v7.3')
end
