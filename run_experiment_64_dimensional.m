%% Set up path
% The following directories need to be subdirectories of adaptive_sc_fem, 
% or independently added to the MATLAB path:
% /sparse-grids-matlab-kit_v-22-02_California
% /ifiss3.6

addpath(genpath(pwd),'-end')

%% Define problem
problem = define_problem('doubleglazing-64');
reference = define_reference('none');

%% Run experiments
names = {'test-interp-jomp'};

for name = names
    params = define_params(name{1});

    adaptive_sc_fem;
    save([name{1} '.mat'],'reference','data_table','fem','problem','params','-v7.3')
end

plot_data(data_table, reference,fem,params,problem)