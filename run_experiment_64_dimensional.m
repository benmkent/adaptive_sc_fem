%% Define problem
problem = define_problem('doubleglazing-64');

%% Run experiments
names = {'test-interp-jomp'};


for name = names
    params = define_params(name{1});

    adaptive_sc_fem;
    save([name{1} '.mat'],'reference','data_table','fem','problem','params','-v7.3')
end

plot_data(data_table, reference,fem,params,problem)