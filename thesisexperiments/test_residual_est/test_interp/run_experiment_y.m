
problem = define_problem('doubleglazing');

% %% Set up reference
% params = define_params('test-interp-ref');
% params.jomp_estimator = 0;
% params.residual_estimator = 0;
% 
% adaptive_sc_fem;
% save(['reference.mat'],'reference','data_table','fem','problem','params')

%% Run experiments
names = {'test-interp'};

problem = define_problem('doubleglazing');

for l_initial = 1:4
    reference = define_reference('test-folder');
    params = define_params('test-interp');
    
    params.jomp_estimator = 0;
    params.residual_estimator = 1;

    params.l_initial = l_initial;

    adaptive_sc_fem;
    save(['interp' num2str(l_initial) '.mat'],'reference','data_table','fem','problem','params')
end