function [error, Z_I_star, precompute] = error_estimate_bk(Z_I_star, Z_I_star_lofi, I_star, I, params, fem, precompute);
fprintf('.... Compute pi_I')
pi_I = compute_pi_I(I_star, I, precompute);
fprintf('...%f\n', pi_I);
fprintf('.... Compute pi_delta, pi_I_delta')
[pi_delta, pi_I_delta, Z_I_star] = compute_pi_delta(Z_I_star, Z_I_star_lofi, I_star, I, params, precompute);
fprintf('...%f...%f\n', pi_delta, pi_I_delta);
pi = pi_I + pi_delta + pi_I_delta;

error.pi = pi;
error.pi_I = pi_I;
error.pi_delta = pi_delta;
error.pi_I_delta= pi_I_delta;
error.pi_I_alpha = {};