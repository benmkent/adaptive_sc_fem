function [flowx,flowy] = gauss_transprt_bk(s,t,xl,yl,wind_fn)
%GAUSS_TRANSPRT evaluates convection field at Gauss point
%   [flowx,flowy] = gauss_transprt(s,t,xl,yl);
%   input
%          s         reference element x coordinate
%          t         reference element y coordinate
%          xl        physical element x vertex coordinates
%          yl        physical element y vertex coordinates
%
%   calls function: specific_wind
%   IFISS function: DJS; 4 March 2005.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage
nel=length(xl(:,1));
nvtx = size(xl,2);
zero_v = zeros(nel,1); xx=zero_v; yy=xx;
[phi_e,dphids,dphidt] = shape(s,t);
for ivtx=1:nvtx
    xx = xx + phi_e(ivtx) * xl(:,ivtx);
    yy = yy + phi_e(ivtx) * yl(:,ivtx);
end
%       [flowx,flowy] = specific_wind(xx,yy,nel);
      flow = wind_fn(xx(:),yy(:),nel);
      flowx = flow(:,1);
      flowy = flow(:,2);
% for ii = 1:nel
%     flow_ii = wind_fn(xx(ii),yy(ii),nel);
%     flowx(ii) = flow_ii(1);
%     flowy(ii) = flow_ii(2);
% end
return
