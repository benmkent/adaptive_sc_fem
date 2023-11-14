function [eta_t_z,eta_t_k_z, eta_t_z_du, eta_t_z_u_w, eta_t_z_eta_w] = compute_temporal_errors(z, fem,tr,trp1)
%COMPUTE_TEMPORAL_ERRORS Computes temporal error estimators
%
% Inputs    z       data structure for a collocation point
%           fem     data structure for FEM approximation
%           tr      time interval start
%           trp1    time interval end
%
% Outputs   eta_z   spatial error estimator
%           eta_z_T localised spatial error estimator on elements

%% Set up problem
theta = 0.5;
xy = fem.xy;
ev = fem.ev;
ebound = fem.ebound;
stiffness = fem.stiffness;
mass = fem.mass;
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
N = z.N;
Q = z.Q;

S = eye(size(Q)); S(bound,:) = [];

fx = @(t) S *f - A(notbound,bound)*bc(t) - Q(notbound,bound)*bc_prime(t);
f_K(bound) = 0;

%% Set up local error problems
% Computing edge lengths/connections
evt = ev;
[eex,tve,els] = tedgegen(xy,evt);
eboundt = ebound;

% Preallocate
n_k = length(z.t_z);
u_w_k_z = zeros(size(U,1),n_k-1);

% Do we set BC to zero?
U(bound,:) = 0;

% Loop over timesteps
for ii_k = 1 : (n_k-1)
    dt_k = DT(ii_k+1);
    dU(:,ii_k) = (U(:,ii_k+1) - U(:,ii_k));

    % Hard code BC fn
    thetaT = theta * time(ii_k+1)+ (1-theta) * time(ii_k);
    bc_fn = @(x1,x2) (1-exp(-10*thetaT))*(1-x2.^4).*(abs(x1-1) < 1e-6);

    % Compute $\Vert \nabla (Ukp1 - Uk) \Vert_{H1}$
    dU_k_z(ii_k) = sqrt(dU(:,ii_k).' * mass * dU(:,ii_k));

    % Compute solution to (u,v)_H^1_0 = ( w \dot dU , v)_H^{-1}
    u_w_k_z(notbound,ii_k) = stiffness(notbound,notbound) \ N(notbound,notbound) * dU(notbound,ii_k);
    u_w_k_z_normed(ii_k) = sqrt(u_w_k_z(:,ii_k).' * stiffness  * u_w_k_z(:,ii_k));

    % Compute a posteriori error estimation for u_w_k_z
    fprintf('Error estimation using 4 bubble functions\n');
    [elerr_nobc,fe,ae] = diffpost_p1_with_p2_for_wind(xy,evt,eex,tve,els,eboundt,u_w_k_z(:,ii_k),dU(:,ii_k),w_fn);
    % zero BC
    bc_fn = @(x1,x2) 0*x1;
    [err_p,elerr(:,ii_k)] = cdpost_p1_bc(ae,fe,elerr_nobc,xy,evt,eboundt,bc_fn);

    % Global error estimate
    errest(ii_k) = norm(elerr,2);

    eta_t_k_z(ii_k) = sqrt(dU_k_z(ii_k)^2 + u_w_k_z_normed(ii_k)^2 + errest(ii_k)^2);
end

%% Adjust 1/12 factor for end points (possible partial integration over time intervals)
if length(time)>2
    t_integral = (1/12).*DT(2:end);
    t_integral(1) = t_integral(1) - (0.25*(tr-time(1)) - 0.5*(DT(2)^-1) * (tr - time(1))^2+ 1/3 * (DT(2))^-2 * (tr - time(1))^3);
    t_integral(end) =(0.25*(trp1-time(end-1)) - 0.5*(DT(end)^-1) * (trp1-time(end-1))^2+ 1/3 * (DT(end))^-2 * (trp1-time(end-1))^3);
else
    t_integral(1) = (0.25*(trp1-time(end-1)) - 0.5*(DT(end)^-1) * (trp1-time(end-1))^2+ 1/3 * (DT(end))^-2 * (trp1-time(end-1))^3)-...
        (0.25*(tr-time(1)) - 0.5*(DT(2)^-1) * (tr - time(1))^2+ 1/3 * (DT(2))^-2 * (tr - time(1))^3);
end

%% Compute total estimator
eta_t_z_du = sqrt(sum(dU_k_z(:).^2 .*t_integral));
eta_t_z_u_w = sqrt(sum(u_w_k_z_normed(:).^2 .*t_integral));
eta_t_z_eta_w = sqrt(sum(errest(:).^2 .*t_integral));

eta_t_z = sqrt(sum(eta_t_k_z(:).^2 .*t_integral));
end