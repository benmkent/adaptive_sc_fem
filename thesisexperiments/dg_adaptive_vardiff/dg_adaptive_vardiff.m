addpath ../../adaptive_sc_fem/;
problem = define_problem('doubleglazing-vardiff');

%% Set up reference
params = define_params('smolyak-ref');
adaptive_sc_fem;
save(['reference.mat'],'reference','data_table','fem','problem','params', '-v7.3')

%% Run experiments
% names = {'smolyak-1','smolyak-2','smolyak-3','smolyak-4','smolyak-5'};
% names = {'smolyak-1','smolyak-2','smolyak-3','smolyak-4','smolyak-5'};
% names = {'smolyak-3','smolyak-4','smolyak-5'};
% names = {'smolyak-1','smolyak-2'};

% problem = define_problem('doubleglazing-206');

    reference = define_reference('test-folder');
    params = define_params('adaptive-base');
%     params.plot = 1;

    adaptive_sc_fem;
    save(['adaptive.mat'],'reference','data_table','fem','problem','params','-v7.3')