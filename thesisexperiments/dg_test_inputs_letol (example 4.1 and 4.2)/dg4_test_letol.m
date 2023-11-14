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

params.letol = 1e-6; % Local error tolerance

adaptive_sc_fem;
save(['l4-jomp-1e6.mat'],'reference','data_table','fem','problem','params','-v7.3')
cd('1em6');
% plot_data(data_table, reference,fem,params,problem)
cd('..');

%%
clearvars, close all
reference = define_reference('test-folder');
params = define_params('adaptive-base');

params.letol = 1e-5; % Local error tolerance

adaptive_sc_fem;
save(['l4-jomp-1e5.mat'],'reference','data_table','fem','problem','params','-v7.3')
cd('1em5');
% plot_data(data_table, reference,fem,params,problem)
cd('..');

%%
clearvars, close all
reference = define_reference('test-folder');
params = define_params('adaptive-base');

params.letol = 1e-4; % Local error tolerance


adaptive_sc_fem;
save(['l4-jomp-1e4.mat'],'reference','data_table','fem','problem','params','-v7.3')
cd('1em4');
% plot_data(data_table, reference,fem,params,problem)
cd('..');

%%
clearvars, close all
reference = define_reference('test-folder');
params = define_params('adaptive-base');

params.letol = 1e-3; % Local error tolerance


adaptive_sc_fem;
save(['l4-jomp-1e3.mat'],'reference','data_table','fem','problem','params','-v7.3')
cd('1em3');
% plot_data(data_table, reference,fem,params,problem)
cd('..');
