function [intres] = cd_intres_p1_with_p2(fem, xy,evt,p1sol, p1timederiv, w_fn, time)
%INTRES_P1_WITH_P2 interior residuals for P1 solution using P2 bubble functions
%
%   [intres] = cd_intres_p1_with_p2(xy,xl_s,yl_s,evt,p1sol)
%
%   input:
%          xy         vertex coordinate vector
%          evt        element mapping matrix
%          p1sol      vertex solution vector
%          p1timederiv      vertex solution vector time derivative
%
%   output:
%          intres     interior residuals
%
% Function(s) called: triangular_gausspoints
%                     tgauss_gradcoeff
%                     tgauss_source
%
% See also INTRES_P1_WITH_P1
%
% Last update: 01/02/2017
% -------------------------------------------------------------------------
%    TIFISS function:
% Copyright (c) 2017 Alex Bespalov, Leonardo Rocchi, modifed BMK 2023

x = xy(:,1);
y = xy(:,2);
nel = length(evt(:,1));

% Recover local coordinates and local solution
% Additionally require local BC function to compute residuals
bc = zeros(size(p1sol));
bc_prime = bc;

bc(fem.bound) = fem.bc_fn(time);
bc_prime(fem.bound) = fem.bc_fn_prime(time);

for ivtx = 1:3
    xl_v(:,ivtx) = x(evt(:,ivtx));
    yl_v(:,ivtx) = y(evt(:,ivtx));
    sl_v(:,ivtx) = p1sol(evt(:,ivtx));
    sldt_v(:,ivtx) = p1timederiv(evt(:,ivtx));
    bc_v(:,ivtx) = bc(evt(:,ivtx));
    bc_prime_v(:,ivtx) = bc_prime(evt(:,ivtx));
end

% Construct 2D gaussian rule over the reference triangle
% nngpt = 7;
% [s,t,wt] = triangular_gausspoints(nngpt);
% USe centroid approximation
nngpt = 1; 
s(1)=0.333333333333333;  t(1)=0.333333333333333;  wt(1)=1;

% Preallocate matrices
intres = zeros(nel,4);
bde = zeros(nel,4,3);
fde = zeros(nel,4);

% Loop over Gauss points
for igpt = 1:nngpt
    sigpt = s(igpt);
    tigpt = t(igpt);
    wght = wt(igpt);

    % Evaluate derivatives
    [jac_v,invjac_v,phi_v,dphidx_v,dphidy_v] = tderiv(sigpt,tigpt,xl_v,yl_v);
    [psi_v,~,~] = tqderiv(sigpt,tigpt,xl_v,yl_v);

    % Gradient of the diffusion coefficients. CURRENTLY FIXED TO ZERO
    [diffx,diffy] = tgauss_gradcoeff(sigpt,tigpt,xl_v,yl_v);

    % Source f. Zero forcing in BMK test problem.
    [rhs_f] = 0*tgauss_source(sigpt,tigpt,xl_v,yl_v);

    %rhs_f = gauss_bc(sigpt, tigpt,xl_v,yl_v, fem,time);

    nel = length(xl_v(:,1));
    zero_v = zeros(nel,1);
    xx = zero_v;
    yy = xx;

    % Wind field at integration pt
    [phi_e,dphids,dphidt] = tshape(sigpt,tigpt);
    for ivtx=1:3
        xx = xx + phi_e(ivtx) * xl_v(:,ivtx);
        yy = yy + phi_e(ivtx) * yl_v(:,ivtx);
    end
    [flow] = w_fn(xx,yy);
    windx = flow(:,1);
    windy = flow(:,2);

    % Loop over the four bubble functions
    for j = 1:4

        % Compute rhs-contribution from the source f
        fde(:,j) = fde(:,j) + wght * rhs_f(:) .* psi_v(:,j+3) .* jac_v(:);

        % Compute div(a*grad)-contribution = grad(a)*grad(u_tau): loop
        % over vertices hat functions
        for i = 1:3
            bde(:,j,i) = bde(:,j,i) + wght * diffx(:) .* (bc_v(:,i) - sl_v(:,i)).* dphidx_v(:,i) .* psi_v(:,j+3);
            bde(:,j,i) = bde(:,j,i) + wght * diffy(:) .* (bc_v(:,i) - sl_v(:,i)).* dphidy_v(:,i) .* psi_v(:,j+3);
            %         end
            % Compute - w dot \grad u
            %         for i = 1:3
            bde(:,j,i) = bde(:,j,i) - wght * windx(:).*(bc_v(:,i) + sl_v(:,i)) .* dphidx_v(:,i) .* psi_v(:,j+3);
            bde(:,j,i) = bde(:,j,i) - wght * windy(:).*(bc_v(:,i) + sl_v(:,i)).* dphidy_v(:,i)  .* psi_v(:,j+3);
            % Compute  - Dt U
            bde(:,j,i) = bde(:,j,i) - wght .* (bc_prime_v(:,i) + sldt_v(:,i)) .* phi_v(:,i) .*jac_v(:) .* psi_v(:,j+3);
        end
        % end vertices hat functions loop
    end
    % end four bubble functions loop
end
% end Gauss points loop

% Assemble interior residuals from rhs-contributions
for i = 1:4
    intres(:,i) = intres(:,i) + fde(:,i);
end

% Multiply div(a*grad)-contribution by Galerkin solution
for j = 1:4
    for k = 1:3
        intres(:,j) = intres(:,j) + bde(:,j,k);
    end
end

end  % end function
