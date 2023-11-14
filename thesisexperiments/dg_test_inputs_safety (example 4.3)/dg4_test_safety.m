%% Define problem
problem = define_problem('doubleglazing');

%% Set up reference for assessing error estimates
params = define_params('smolyak-ref');
reference = define_reference('none');

adaptive_sc_fem;
save(['reference.mat'],'reference','data_table','fem','problem','params', '-v7.3')

%% Run experiments

clearvars, close all
reference = define_reference('test-folder');
params = define_params('adaptive-base');

params.k_interp = 2; % Safety factor

adaptive_sc_fem;
save(['l4-jomp-k2.mat'],'reference','data_table','fem','problem','params','-v7.3')

%%
clearvars, close all
reference = define_reference('test-folder');
params = define_params('adaptive-base');

params.k_interp = 10^2; % Safety factor

adaptive_sc_fem;
save(['l4-jomp-k100.mat'],'reference','data_table','fem','problem','params','-v7.3')


