function [error, Z_I_star] = error_estimate_bk(Z_I_star, Z_I_star_lofi, I_star, I, params, precompute)
%ERROR_ESTIMATE_BK  Compute hierarchical error estimates
%
% Inputs        Z_I_star        data structure of collocation points
%               Z_I_star_lofi   lofi data structure of collocation points
%               I_star          Enhanced sparse grid
%               I               Sparse grid
%               params          Approximation parameters
%               precompute      Precomputed data
% Outputs       error           Structure of computed error estimates
%               Z_I_star        Updated data structure

%% Compute estimators
fprintf('.... Compute pi_I')
pi_I = compute_pi_I(I_star, I, precompute);
fprintf('...%f\n', pi_I);
fprintf('.... Compute pi_delta, pi_I_delta')
[pi_delta, pi_I_delta, Z_I_star] = compute_pi_delta(Z_I_star, Z_I_star_lofi, I_star, I, params, precompute);
fprintf('...%f...%f\n', pi_delta, pi_I_delta);
pi = pi_I + pi_delta + pi_I_delta;

%% Return as structure
error.pi = pi;
error.pi_I = pi_I;
error.pi_delta = pi_delta;
error.pi_I_delta= pi_I_delta;
error.pi_I_alpha = {};