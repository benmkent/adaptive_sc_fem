function [diffx,diffy] = tgauss_adiff(s,t,xl,yl,diff)
%TGAUSS_ADIFF  evaluates permeability field at triangle Gauss point
%   [diffx,diffy] = tgauss_adiff(s,t,xl,yl);
%   input
%          s         reference element x coordinate   
%          t         reference element y coordinate
%          xl        physical element x vertex coordinates 
%          yl        physical element y vertex coordinates  
%
%   calls function: specific_adiff
%    PIFISS function: DJS; 31 January 2007. 
% Copyright (c) 2007 C.E. Powell, D.J. Silvester
      nel=length(xl(:,1));
      zero_v = zeros(nel,1); xx=zero_v; yy=xx;
      [xi,dxids,dxidt] = tshape(s,t);
      for ivtx=1:3
      xx = xx + xi(ivtx) * xl(:,ivtx);
      yy = yy + xi(ivtx) * yl(:,ivtx);
	  end
      [diff_xxyy] = diff(xx,yy);  
      diffx=diff_xxyy(:,1); 
      diffy=diff_xxyy(:,4); 
      return