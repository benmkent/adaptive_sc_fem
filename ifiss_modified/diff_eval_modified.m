function [diff_eval] = diff_eval_bk(s,t,xl,yl,diff_fn)
%diff_eval_bk evaluates diffusion field at Gauss point
%   [diff_eval] = diff_eval_bk(s,t,xl,yl,diff_fn);
%   input
%          s         reference element x coordinate
%          t         reference element y coordinate
%          xl        physical element x vertex coordinates
%          yl        physical element y vertex coordinates
%          diff_fn   Diffusion field function
%
%   calls function: specific_wind
%   IFISS function: DJS; 4 March 2005.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage
nel=length(xl(:,1));
zero_v = zeros(nel,1); xx=zero_v; yy=xx;
[phi_e,dphids,dphidt] = shape(s,t);
for ivtx=1:4
    xx = xx + phi_e(ivtx) * xl(:,ivtx);
    yy = yy + phi_e(ivtx) * yl(:,ivtx);
end

diff_eval = diff_fn(xx(:),yy(:));
return
