function [eta_z,eta_z_T] = compute_element_spatial_errors(z, fem,tr,trp1)
theta = 0.5;
epsilon = fem.lambda;
r = 0;

hx = fem.hx;
hy = fem.hy;
eex = fem.eex;
xy = fem.xy;
ev = fem.ev;
ebound = fem.ebound;
tve = fem.tve;

bound = fem.bound;
notbound = fem.notbound;

f = z.f;
U = z.u_z;
DT = z.dt_z;
time = z.t_z;
a_fn = z.a_fn;
w_fn = z.wind_fn;
bc = z.bc;
bc_prime = z.bc_prime;
A = (z.D + z.N);
Q = z.Q;

S = eye(size(Q)); S(bound,:) = [];

fx = @(t) S *f - A(notbound,bound)*bc(t) - Q(notbound,bound)*bc_prime(t);
f_K(bound) = 0;

%% Set up local error problems
% Computing edge lenghts/connections
evt = ev;
[eex,tve,els] = tedgegen(xy,evt);
eboundt = ebound;

n_k = length(z.t_z);
for ii_k = 1 : (n_k-1)
    dt_k = DT(ii_k+1);
    dtU(:,ii_k) = dt_k^(-1) * (U(:,ii_k+1) - U(:,ii_k));
    thetaU(:,ii_k) = theta*U(:,ii_k+1)+ (1-theta)*U(:,ii_k);
    
    % Hard code BC fn
    thetaT = theta * time(ii_k+1)+ (1-theta) * time(ii_k);
    bc_fn = @(x1,x2) (1-exp(-10*thetaT))*(1-x2.^4).*(abs(x1-1) < 1e-6);

    % Compute a posteriori error estimation
    fprintf('Error estimation using 4 quadratic bubble functions\n');
    [elerr_nobc,fe,ae] = cdpost_p1_with_p2(xy,evt,eex,tve,els,eboundt,thetaU(:,ii_k),dtU(:,ii_k),w_fn,a_fn);
    [err_p,elerr(:,ii_k)] = cdpost_p1_bc(ae,fe,elerr_nobc,xy,evt,eboundt,bc_fn);

    % Global error estimate
    errest(ii_k) = norm(elerr,2);
end

% Correct first and final integrals (over only part of each timestep)
if length(time) > 2
    DT(2) = time(2)-tr;
    DT(end) = trp1 - time(end-1);
else
    DT(2) = trp1-tr;
end
eta_z = sqrt(sum(errest(:).^2 .* DT(2:end)));


for ii_T = 1:size(elerr,1)
    eta_z_T(ii_T) = sqrt(sum((elerr(ii_T,:).').^2 .* DT(2:end)));
end
end