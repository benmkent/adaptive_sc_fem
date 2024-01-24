function [wind_fn, lambda, Vl, w_max, w_mean, xp, yp] = djs_wind_perturbation(n_rv, cov, x_level)
%DJS_WIND_PERTURBATION Generate n_rv wind fields through a streamd
%function and fixed covariance function
%   
% Inputs
%   n_rv        number of parameters
%   cov         coefficient of variation
%   x_level     level number for uniform, rectangular spatial mesh for
%               computing stream function, 2^l x^l elements.
% Outputs
%   wind_fn     cell array of wind field functions
%   lambda      computed eigenvalues for covariance function
%   Vl          computed eigenvectors scaled by evalues
%   w_max       approximation of maximum wind magnitude on domain
%   w_mean      approximation of mean wind magnitude on domain
%   xp          x grid for computing evalues and evectors
%   yp          y grid for computing evalues and evectors


% Ben (Catherine, for information)
% As discussed yesterday I am attaching a bunch of m-files
% Note that these are designed to work with Q2 approximation not Q1
%
% ---------------------------------------------
% Start with linstab_startdata.m lines 19 to 51
% line 20 ...
% xyp are the Q1 coordinates of the current grid defined by xy and mv
% These are the vertices coordinates
%[xyp,mp,map] = q2q1map(xy,mv);
% x_level=7;
[ ~, xyp, boundp] = rsquare_domain_modified_fn(x_level,1);
np = length(xyp);

% lines 21 to 38
% sets up discrete covariance matrix Cvar
% (setting delta=1 generates rapidly decaying eigenvalues ...)

% Set up covariance operator
dx = abs( xyp(:,1)*ones(1,np) - ones(np,1)*xyp(:,1)' );
dy = abs( xyp(:,2)*ones(1,np) - ones(np,1)*xyp(:,2)' );

boundary_fix = (1-abs(xyp(:,1)).^(2)).*(1-abs(xyp(:,2)).^(2));
boundary_fix = boundary_fix * boundary_fix';

%corrL = 2.0;    %% reference length L for the domain
%dist = sqrt( dx.^2 + dy.^2 );
%delta = 1;
%Cvar = exp(-(dist./corrL).^(1+delta));

delta = 1; %problem.delta; % delta = 1;
corrx = 1; %2; % problem.corrx; %         corrx = 2;
corry = 1; %2; %problem.corry; %corry = 2;
dist = sqrt((dx/corrx).^2 + (dy/corry).^2);
clear dx dy;
Cvar = exp(-dist.^(1+delta));
clear dist;
Cvar_fix = boundary_fix.* Cvar;
%Cvar_fix = Cvar;
clear Cvar boundary_fix;

% Cvar = min((xyp(:,1)+1)*0.5,2) - prod((xyp(:,1)+1)*0.5,2);

%dist = dx/corrx + dy/corry;
%Cvar = exp(-dist);

% line 43
% cov is a multiplier ("coefficient of variation")

% cov   = 1e-2; %problem.cov;
Cvar_fix = cov * Cvar_fix;

% lines 47 to 49
% generates kl biggest eigenvalues and corresponding eigenvectors
% (captures 95% of total variance)

% Compute Q
h = 2 / 2^(x_level);
Q = spdiags(h^2*ones(size(Cvar_fix,1),1),0,size(Cvar_fix,1),size(Cvar_fix,1));
Q(boundp,boundp) = spdiags(h^2/2 * ones(length(boundp),1),0,length(boundp),length(boundp));
corners = find((abs(xyp(:,1)) == 1) & (abs(xyp(:,2)) == 1));
Q(corners,corners) = spdiags(h^2/4 * ones(length(corners),1),0,length(corners),length(corners));

fprintf('Computing eigenvalues of covariance operator ... ');
CQ = Cvar_fix * Q;
% CQ = Cvar_fix;
Qhalf = spdiags(sqrt(diag(Q)),0);
S = Qhalf * Cvar_fix * Qhalf;
[Z,D] = eigs(S, n_rv);
% [V,D] = eigs(CQ,n_rv,'largestabs','IsFunctionSymmetric','true');
% D = real(D);
[lambda,index] = sort(diag(D),'descend');
% V = V(:,index);
Z = Z(:,index);
V = Qhalf \ Z;
normV = sqrt(diag(Q).' * V.^2);
V = V * diag(1./normV);
fprintf('done.\n');
% Use enough terms to capture 95% of fluctuation
kl = np;
lt = sum(lambda);
% lp = lambda(np);
%         while lp/lt < .05,
%             kl = kl-1;
%             lp = lp + lambda(kl);
%         end
kl = n_rv;
% ratio = sum(lambda(1:kl))/sum(lambda);

% lines 77 to 81
% saves data for later use in truncated "KL"  expansion
% (note the scaling of the eigenvector in line 80
% by the square root of the associated eigenvalue ..)
%         Vl = V(:,1:kl)*diag(sqrt(lambda(1:kl)));
Vl = V(:,1:kl)*diag(sqrt(lambda(1:kl)));

% ---------------------------------------------
% Look at linstab_alleig_for_spinterp.m next
%
% lines 51 to 54
% generate the scalar perturbation field dpsi
% associated the multidimensional collocation point xi
% construct the perturbed veocity field and convection matrix
U = 1;      % normalization of velocity via max(inflow) / max(cavity b.c.)
boundvec = ones(size(xyp(:,1)));
%         boundvec(boundp) = 0;
%         dpsi = @(xi) U*cov*(Vl*xi(:)) .* boundvec;

%         psi0 = xyp(:,1).^2 .* (1-xyp(:,2).^2) + xyp(:,2).^2;

n_x = sqrt(size(xyp,1));
xp = reshape(xyp(:,1), [n_x,n_x]);
yp = reshape(xyp(:,2), [n_x,n_x]);

Vl_xyp = reshape(Vl, [n_x,n_x,n_rv]);

%         wx_xyp = diff(Vl_xyp,2,1)./(diff(xp,2,1));
dphidx_xyp = (Vl_xyp(3:end,:,:)-Vl_xyp(1:(end-2),:,:))./(xp(3:end,:)-xp(1:(end-2),:));
% Extend on to x=-1,1 boundary
dphidx_xyp_extended = zeros(size(Vl_xyp));
dphidx_xyp_extended(2:end-1,:,:) = dphidx_xyp;
%         wx_xyp = wx_xyp(:,2:end-1,:);
%         xp = xp(2:(end-1),2:(end-1),:);
%         wy_xyp = diff(Vl_xyp,2,1)./diff(yp,2,2);
dphidy_xyp = (Vl_xyp(:,3:end,:)-Vl_xyp(:,1:(end-2),:))./(yp(:,3:end)-yp(:,1:(end-2)));
%Extend onto y=-1,1 boundary
dphidy_xyp_extended = zeros(size(Vl_xyp));
dphidy_xyp_extended(:,2:end-1,:) = dphidy_xyp;

%%
wx = dphidy_xyp_extended;
wy = -dphidx_xyp_extended;
w_mag = sqrt(wx.^2 + wy.^2);
w_max=max(squeeze(max(w_mag,[],1)),[],1);
w_mean=mean(squeeze(mean(w_mag,1)),1);

for ii = 1:n_rv
    wind_fn{ii} = @(x1,x2, nel) [griddata(xp,yp,wx(:,:,ii), x1,x2,'cubic'),griddata(xp,yp,wy(:,:,ii), x1,x2,'cubic')];
end

%         wy_xyp = wy_xyp(2:end-1,:,:);
%         yp = yp(2:(end-1),2:(end-1),:);


% lines 62 to 63
% plots the perturbation field
% Note that dpsi will be interpreted as a discrete stream function field
% in what follows
%
% ---------------------------------------------
% Look at linstab_neig_for_spinterp.m next
%
% lines 64 to 68
% construct the perturbation (note multiplication by cov  in line 68)

% line 74
% generate the discrete convection matrix Nxi associated with the
% perturbation (should be skew-symmetric by construction, see below)
%         psi = @(xi) psi0 + dpsi(xi);