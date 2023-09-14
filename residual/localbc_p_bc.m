function [ae,fe] = localbc_p(ae,fe,edges,xl,yl,bc_fn)
%LOCALBC_P  imposes Dirichlet BC for Poisson error estimator
%
%   [ae,fe] = localbc_p(ae,fe,edges,xl,yl)
%
%   input
%               ae     Poisson problem matrix
%               fe     rhs vector
%            edges     boundary edge vector 
%           xl, yl     vertex coordinates  
%   output 
%               ae     Poisson problem matrix
%               fe     rhs vector
%
% Function(s) called:  specific_bc
%
%   TIFISS function: LR; 05 October 2017.
% Copyright (c) 2017 A. Bespalov, L. Rocchi

  nbd = length(edges);  % number of boundary edges (either 1 or 2) 

% Loop over boundary edges
  for bd = 1:nbd
      ek = edges(bd);   % edge's number (1,2, or 3)
      % Recover boundary edge coordinates
      xbd(1) = xl(ek+1);    xbd(3) = xl(ek+2);    xbd(2) = 0.5*(xbd(1) + xbd(3));
      ybd(1) = yl(ek+1);    ybd(3) = yl(ek+2);    ybd(2) = 0.5*(ybd(1) + ybd(3));      
      % Boundary error
      [bc] = bc_fn(xbd,ybd);
      error = bc(2) - 0.5*(bc(1) + bc(3));
      % Impose boundary condition without modifying the other equations (DJS/DK mod)
      % fe = fe - error*ae(:,ek);
      ae(:,ek)  = 0.0;
      ae(ek,:)  = 0.0;
      ae(ek,ek) = 1.0; 
      fe(ek)    = error; 
  end
  
end % end function