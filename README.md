This is the code corresponding to the paper [Efficient Adaptive Stochastic Collocation Strategies for Advection-Diffusion Problems with Uncertain Inputs
](https://link.springer.com/article/10.1007/s10915-023-02247-w) ([Preprint](https://arxiv.org/abs/2210.03389)).

This repository also contains code for the experiments contained within Benjamin M Kent's thesis at the University of Manchester.

# Contents
- [Required MATLAB packages](#required-matlab-packages)
- [Running the examples for the paper](#running-the-examples-for-the-paper)
- [Running the thesis experiments](#running-the-thesis-experiments)

# Required MATLAB packages
To use the package ensure the path contains:
- [IFISS 3.6](https://personalpages.manchester.ac.uk/staff/david.silvester/ifiss/)
- [TIFISS 1.2](https://personalpages.manchester.ac.uk/staff/david.silvester/ifiss/tifiss.html)
- [Sparse Grids MATLAB Kit 22-2 ("California")](https://sites.google.com/view/sparse-grids-kit)

This can be achieved by downloading the appropriate folders are ensuring the directories are placed within the ```adaptive_sc_fem/``` directory. 
These specific versions are required. Function names in the more recent release of the Sparse Grids MATLAB Kit have not been updated and tested in this code yet.

# Running the examples for the paper
There are two examples included corresponding to the examples in [Efficient Adaptive Stochastic Collocation Strategies for Advection-Diffusion Problems with Uncertain Inputs
](https://link.springer.com/article/10.1007/s10915-023-02247-w).
These are the scripts
```matlab
run_experiment_4_dim.m
```
and
```matlab
run_experiment_64_dimensional.m
```
To save on computational time it is possible to run ```run_experiment_4_dim_noreference.m``` and only construct an approximation.

The output results from the four parameter problem are *slightly* different when compared to those plotted in the paper due to small updates to the space--time approximation implementation.

The four parameter example is annotated below.

Firstly the path is set up.
Folders containing ```ifiss3.6```, ``tifiss1.2``` and ```sparse-grids-matlab-kit``` must be either be in the current directory, or added to the path by the user.
```matlab
%% Set up path
% The following directories need to be subdirectories of adaptive_sc_fem, 
% or independently added to the MATLAB path:
% /sparse-grids-matlab-kit_v-22-02_California
% /ifiss3.6

addpath(genpath(pwd),'-end')
```

The test problem is then defined.
```matlab
%% Define problem
problem = define_problem('doubleglazing');
```
A reference solution is then set up. Approximation errors will be computed with respect to this reference solution. The main loop is run and produces a high fidelity reference approximation.
```matlab
%% Set up reference for assessing error estimates
params = define_params('l4-jomp');
reference = define_reference('none');

params.l_initial = 5; % Set reference Smolyak sparse grid level (polynomials including TD 5).
params.letol = 1e-7; % Set reference local error tolerance
params.reference = 1; % Set as reference approximation
params.k_interp = inf; % Set parametric refinement threshold to inf (no refinement)
adaptive_sc_fem;
save(['reference.mat'],'reference','data_table','fem','problem','params', '-v7.3')
```

The approximation is then constructed.
The reference approximation is specified to be precomputed and saved within the current folder (the ```test-folder``` option).
The parameters are set to the predefined set ```l4-jomp```.
```matlab
%% Run experiments

reference = define_reference('test-folder');
params = define_params('l4-jomp');
```
Approximation parameters can then be varied within the ```params``` structure.
```
params.adapt_interp = 0; % Interpolate for new collocation points
params.k_interp = 10; % c_{tol} safety factor
params.marking_factor = 0.1; % Marking factor
params.letol = 1e-5; % Local error tolerance
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
For example, the resulting error estimate, error and tolerance are output as the following figure.

![Error, error estimate and tolerance for four parameter test problem](https://github.com/benmkent/adaptive_sc_fem/assets/52756911/9fad03e2-e509-4dc6-ade1-5bdf35196259)

# Running the thesis experiments
The thesis experiments can be found within subfolders the folder ```thesisexperiments```. The folder names are representative of the experiments within.
The path must be set up as described [above](#running-the-examples-for-the-paper).
Examples in the thesis without dedicated folders are either easily recreated or the plotted data is derived from data generated in other experiments.
The output ```*.mat``` files will generally need postprocessing with the ```subroutine/plot_data.m``` function.