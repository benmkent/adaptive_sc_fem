function E = set_threshold_E(error, params)
%SET_THRESHOLD_E    Defines the interpolation error tolerance E
%
%   Inputs  error   error structure
%           params  approximation parameters
%   Outputs E       interpolation error tolerance

pi_delta = error.pi_delta;
pi_I_delta = error.pi_I_delta;
k_interp = params.k_interp;

E = max([pi_delta,pi_I_delta])*k_interp;