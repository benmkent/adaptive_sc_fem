function E = set_threshold_E(error, params);
pi_delta = error.pi_delta;
pi_I_delta = error.pi_I_delta;
k_interp = params.k_interp;

E = max([pi_delta,pi_I_delta])*k_interp;