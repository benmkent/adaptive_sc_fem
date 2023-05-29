function [n, xl_v, yl_v,flowx, flowy] = femq1_conv(xy,ev, wind_fn)
%FEMQ1_CD vectorized bilinear coefficient matrix generator
%    [N] = femq1_cd(xy,ev);
%   input
%          xy       vertex coordinate vector
%          ev       element mapping matrix
%   output
%          N        convection matrix
%   Natural boundary conditions apply. Dirichlet conditions
%   must be explicitly enforced by calling function nonzerobc.
%   IFISS function: DJS; 5 March 2005.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage


x=xy(:,1); y=xy(:,2);
nvtx=length(x);
nel=length(ev(:,1));
lx=max(x)-min(x); ly=max(y)-min(y);
hx=max(diff(x)); hy=max(diff(y));
fprintf('setting up Q1 convection matrix...\n')
%
% initialise global matrices
n = sparse(nvtx,nvtx);
%
% % set up 2x2 Gauss points
%      ngpt = 4;
%       gpt=1.0e0/sqrt(3.0e0);
%       s(1) = -gpt;  t(1) = -gpt;
%       s(2) =  gpt;  t(2) = -gpt;
%       s(3) =  gpt;  t(3) =  gpt;
%       s(4) = -gpt;  t(4) =  gpt;

% set up 3x3 Gauss points
gpt=sqrt(0.6); ngpt = 9;
s(1) = -gpt; t(1) = -gpt; wt(1)=25/81;
s(2) =  gpt; t(2) = -gpt; wt(2)=25/81;
s(3) =  gpt; t(3) =  gpt; wt(3)=25/81;
s(4) = -gpt; t(4) =  gpt; wt(4)=25/81;
s(5) =  0.0; t(5) = -gpt; wt(5)=40/81;
s(6) =  gpt; t(6) =  0.0; wt(6)=40/81;
s(7) =  0.0; t(7) =  gpt; wt(7)=40/81;
s(8) = -gpt; t(8) =  0.0; wt(8)=40/81;
s(9) =  0.0; t(9) =  0.0; wt(9)=64/81;

%
% inner loop over elements
for ivtx = 1:4
    xl_v(:,ivtx) = x(ev(:,ivtx));
    yl_v(:,ivtx) = y(ev(:,ivtx));
end
ne = zeros(nel,4,4);
% loop over Gauss points
for igpt = 1:ngpt
    sigpt=s(igpt);
    tigpt=t(igpt);
    wigpt = wt(igpt);
    %  evaluate derivatives etc
    [jac,invjac,phi,dphidx,dphidy] = deriv(sigpt,tigpt,xl_v,yl_v);
    %          [flowx,flowy] = gauss_transprt_modified(sigpt,tigpt,xl_v,yl_v,wind_fn);

    nel=length(xl_v(:,1));
    nevtx = size(xl_v,2);
    zero_v = zeros(nel,1); xx=zero_v; yy=xx;
    [phi_e,dphids,dphidt] = shape(sigpt,tigpt);
    for ivtx=1:nevtx
        xx = xx + phi_e(ivtx) * xl_v(:,ivtx);
        yy = yy + phi_e(ivtx) * yl_v(:,ivtx);
    end
    %       [flowx,flowy] = specific_wind(xx,yy,nel);
    flow = wind_fn(xx(:),yy(:),nel);
    flowx = flow(:,1);
    flowy = flow(:,2);

    for j = 1:4
        for i = 1:4
            ne(:,i,j) = ne(:,i,j) + wigpt*flowx(:) .* phi(:,i) .* dphidx(:,j);
            ne(:,i,j) = ne(:,i,j) + wigpt*flowy(:) .* phi(:,i) .* dphidy(:,j);
        end
    end
    % end of Gauss point loop
end
%
% perform assembly of global matrix  and source vector
for krow=1:4
    nrow=ev(:,krow);
    for kcol=1:4
        ncol=ev(:,kcol);
        n = n + sparse(nrow,ncol,ne(:,krow,kcol),nvtx,nvtx);
    end
end

n = 0.5*(n - n.');
return
