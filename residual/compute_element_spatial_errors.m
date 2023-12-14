function [eta_z,eta_z_T] = compute_element_spatial_errors(z, fem,tr,trp1)
%COMPUTE_ELEMENT_SPATIAL_ERRORS Compute thte spatial error estimators
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
% Computing edge lengths/connections
evt = ev;
[eex,tve,els] = tedgegen(xy,evt);
eboundt = ebound;

n_k = length(z.t_z);
% We must now iterate over each time interval $J_{r,k}^{z}$
for ii_k = 1 : (n_k-1)
    % Set up
    dt_k = DT(ii_k+1);
    dtU(:,ii_k) = dt_k^(-1) * (U(:,ii_k+1) - U(:,ii_k));
    thetaU(:,ii_k) = theta*U(:,ii_k+1)+ (1-theta)*U(:,ii_k);
    
    % Hard code BC fn - should use infro from problem or FEM
    thetaT = theta * time(ii_k+1)+ (1-theta) * time(ii_k);
    bc_fn = @(x1,x2) (1-exp(-10*thetaT))*(1-x2.^4).*(abs(x1-1) < 1e-6);

    % Compute a posteriori error estimates solving local error problem
    % $(u,v) = (f,v) - (a\nabla^2 u, v) - (w\nablau,v) -(u_{k+1} - u_k/
    % dt_k,v) - 0.5[[a \nabla u]]$.
    fprintf('Error estimation using 4 quadratic bubble functions\n');
    [elerr_nobc,fe,ae] = cdpost_p1_with_p2(fem,xy,evt,eex,tve,els,eboundt,thetaU(:,ii_k),dtU(:,ii_k),w_fn,a_fn,thetaT);
    % Use TIFISS BC on Dirichlet boundary -- CHECK DO WE ACTUALLY WANT TO
    % DO THIS? This should probably be replaced using the data estimator.
%     [err_p,elerr(:,ii_k)] = cdpost_p1_bc(ae,fe,elerr_nobc,xy,evt,eboundt,bc_fn);
    elerr(:,ii_k) = elerr_nobc; % DO NOT APPLY THE BOUNDARY CORRECTION (THE ESTIMATOR HAS ZERO ERROR AT BOUNDARY)
    
    % Global error estimate
    errest(ii_k) = norm(elerr,2);
end

%% Correct first and final integrals (may only be part of each timestep)
if length(time) > 2
    DT(2) = time(2)-tr;
    DT(end) = trp1 - time(end-1);
else
    DT(2) = trp1-tr;
end
%% Compute final estimates on interval J_{r}
eta_z = sqrt(sum(errest(:).^2 .* DT(2:end)));
for ii_T = 1:size(elerr,1)
    eta_z_T(ii_T) = sqrt(sum((elerr(ii_T,:).').^2 .* DT(2:end)));
end
end