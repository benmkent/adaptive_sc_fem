function [n, xl_v, yl_v,flowx, flowy] = femp1_conv(xy,ev, wind_fn)
%FEMP1_CONV vectorized linear coefficient matrix generator for P1
%convection terms
%
% This function is based upon FEMQ1_CD vectorized bilinear coefficient 
% matrix generator.
%
% BMK 2022
%
%    [N] = femq1_cd(xy,ev);
%   input
%          xy       vertex coordinate vector
%          ev       element mapping matrix
%          wind_fn  advection field to generate matrix for
%   output
%          N        convection matrix
%   Natural boundary conditions apply. Dirichlet conditions
%   must be explicitly enforced by calling function nonzerobc.
%   IFISS function: DJS; 5 March 2005.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage


x=xy(:,1); y=xy(:,2);
nvtx=length(x);
nel=length(ev(:,1));
% lx=max(x)-min(x); ly=max(y)-min(y);
% hx=max(diff(x)); hy=max(diff(y));
fprintf('setting up P1 convection matrix...\n')
%
% initialise global matrices
n = sparse(nvtx,nvtx);

%------------------------------------------------
% 3 point Gauss rule integration
nngpt=7; [s,t,wt]=triangular_gausspoints(nngpt);


% inner loop over elements
for ivtx = 1:3
    xl_v(:,ivtx) = x(ev(:,ivtx));
    yl_v(:,ivtx) = y(ev(:,ivtx));
end

ne = zeros(nel,3,3);

%
% loop over Gauss points
for igpt = 1:nngpt
    sigpt=s(igpt);
    tigpt=t(igpt);
    wtigpt=wt(igpt);
    %  evaluate derivatives etc
    [jac,invjac,phi,dphidx,dphidy] = tderiv(sigpt,tigpt,xl_v,yl_v);
%     [flowx,flowy] = gauss_transprt_modified(sigpt,tigpt,xl_v,yl_v,wind_fn);
    nel=length(xl_v(:,1));
    nevtx = size(xl_v,2);
    zero_v = zeros(nel,1); xx=zero_v; yy=xx;
    [phi_e,dphids,dphidt] = tshape(sigpt,tigpt);
    for ivtx=1:nevtx
        xx = xx + phi_e(ivtx) * xl_v(:,ivtx);
        yy = yy + phi_e(ivtx) * yl_v(:,ivtx);
    end
    %       [flowx,flowy] = specific_wind(xx,yy,nel);
   flow = wind_fn(xx(:),yy(:),nel);
   flowx = flow(:,1);
   flowy = flow(:,2);

    for j = 1:3
        for i = 1:3
            ne(:,i,j) = ne(:,i,j) +  wtigpt * flowx(:) .* phi(:,i) .* dphidx(:,j);
            ne(:,i,j) = ne(:,i,j) + wtigpt * flowy(:) .* phi(:,i) .* dphidy(:,j);
        end
    end
    % end of Gauss point loop
end
%
% perform assembly of global matrix  and source vector
for krow=1:3
    nrow=ev(:,krow);
    for kcol=1:3
        ncol=ev(:,kcol);
        n= n + sparse(nrow,ncol,ne(:,krow,kcol),nvtx,nvtx);
    end
end
return
