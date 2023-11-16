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

params.marking_factor = 0.5; % Local error tolerance

adaptive_sc_fem;
save(['l4-jomp-t05.mat'],'reference','data_table','fem','problem','params','-v7.3')
cd('theta05');
process_output_data(data_table, reference,fem,params,problem)
cd('..');

%%
clearvars, close all
reference = define_reference('test-folder');
params = define_params('adaptive-base');

params.marking_factor = 1e-2; % Local error tolerance

adaptive_sc_fem;
save(['l4-jomp-t001.mat'],'reference','data_table','fem','problem','params','-v7.3')
cd('theta001');
process_output_data(data_table, reference,fem,params,problem)
cd('..');

%%
clearvars, close all
reference = define_reference('test-folder');
params = define_params('adaptive-base');

adaptive_sc_fem;
save(['l4-jomp-t01.mat'],'reference','data_table','fem','problem','params','-v7.3')
cd('theta01');
process_output_data(data_table, reference,fem,params,problem)
cd('..');
