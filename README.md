This is the code corresponding to the preprint [Efficient Adaptive Stochastic Collocation Strategies for Advection-Diffusion Problems with Uncertain Inputs
](https://arxiv.org/abs/2210.03389).

# Required MATLAB packages
To use the pacakge ensure the path contains:
- [IFISS 3.6](https://personalpages.manchester.ac.uk/staff/david.silvester/ifiss/)
- [Sparse Grids MATLAB Kit](https://sites.google.com/view/sparse-grids-kit)

# Running the examples
There are two examples included corresponding to the examples in [Efficient Adaptive Stochastic Collocation Strategies for Advection-Diffusion Problems with Uncertain Inputs
](https://arxiv.org/abs/2210.03389).
These are the scripts
```matlab
run_experiment_4_dim.m
```
and
```matlab
run_experiment_64_dimensional.m
```
The four parameter example is annotated below.

First the test problem is defined.
```matlab
%% Define problem
problem = define_problem('doubleglazing');
```
A reference solution is then set up. Approximation errors will be computed with respect to this reference solution. A prompt stating
```Load file reference.mat :``` will appear. Input ```0``` to continue (this creates a reference solution instead of loading a precomputed one).
```matlab
%% Set up reference
params = define_params('l4-jomp');
params.l_initial = 5;
params.letol = 1e-7;
adaptive_sc_fem;
save(['reference.mat'],'reference','data_table','fem','problem','params', '-v7.3')
```

The approximation is then constructed.
The reference approximation is defined to be precomputed and saved within the current folder.
The parameters are set to the predefined set ```l4-jomp```.
```matlab
%% Run experiments

reference = define_reference('test-folder');
params = define_params('l4-jomp');
```
Approximation parameters can then be varied within the ```params``` structure.
```
params.adapt_interp = 1; % Interpolate for new collocation points
params.simplified_estimator = 1; % Alternative error estimator (all colloc pts use the same)
params.k_interp = 10; % c_{tol} safety factor
params.marking_factor = 0.1; % Marking factor
params.letol = 1e-4; % Local error tolerance
```
The approximation is constructed and the outputs saved to ```l4-jomp.mat```.
```matlab
adaptive_sc_fem;
save(['l4-jomp.mat'],'reference','data_table','fem','problem','params','-v7.3')
```
The results are post processed into MATLAB figures and CSV files.
```matlab
plot_data(data_table, reference,fem,params,problem)
```
