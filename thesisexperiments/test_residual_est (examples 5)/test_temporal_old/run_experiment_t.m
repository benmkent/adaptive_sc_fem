addpath ../../adaptive_sc_fem/;

problem = define_problem('doubleglazing-no-rv');

%% Set up reference 
params = define_params('test-spatial-ref');
params.grid_param = 7;
params.letol=1e-7;

adaptive_sc_fem;
save(['reference.mat'],'reference','data_table','fem','problem','params')

%% Name experiments
names = {'test-temporal-le1e0',...
    'test-temporal-le1e-1',...
    'test-temporal-le1e-2',...
    'test-temporal-le1e-3',...
    'test-temporal-le1e-4',...
    'test-temporal-le1e-5',...
    };
%,...
 %   'test-test-spatial-l6'};

 %% Run experiments
for name = names
reference = define_reference('test-folder');
params = define_params(name{1});
params.residual_estimator = 1;
params.grid_param = 6;

adaptive_sc_fem;
save([name{1} '.mat'],'reference','data_table','fem','problem','params')
end