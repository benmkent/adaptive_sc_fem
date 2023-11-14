problem = define_problem('doubleglazing-no-rv');
reference = define_reference();
params = define_params('test-spatial-l3');

params.plot = 1;
params.grid_param=3;

[data_table,fem,problem,params,reference] = adaptive_sc_fem(problem,params,reference);