function [edgeres] = cd_edgeres_p1_with_p2(xy,evt,eboundt,p1sol,eex,tve,els,diff)
%EDGERES_P1_WITH_P2 edge residuals for P1 solution using P2 bubble functions
%
%   [edgeres] = edgeres_p1_with_p2(xy,evt,eboundt,p1sol,eex,tve,els)
%
%   input:
%               xy    vertex coordinate vector
%              evt    element mapping matrix
%          eboundt    element boundary mapping matrix
%            p1sol    vertex solution vector
%              eex    element connectivity array
%              tve    edge location array
%              els    elementwise edge lengths
%
%   output:
%          edgeres    edge residuals
%
% Function(s) called: gausspoints_oned
%                     p1fluxjmps 
%                     reorder_s
%                     vtqderiv   
%
% See also EDGERES_P1_WITH_P1
%
% Last update: 01/02/2017
% -------------------------------------------------------------------------
%    TIFISS function:
% Copyright (c) 2017 Alex Bespalov, Leonardo Rocchi

  x = xy(:,1);
  y = xy(:,2);
  nel = length(evt(:,1));
   
% Recover local coordinates and solution   
  for ivtx = 1:3 
      xl_v(:,ivtx) = x(evt(:,ivtx));
      yl_v(:,ivtx) = y(evt(:,ivtx));
      sl_v(:,ivtx) = p1sol(evt(:,ivtx));
  end

% Construct the 1D integration rule
  ngpt = 7;
  [oneg,onew] = gausspoints_oned(ngpt);

% Preallocate matrix
  edgeres = zeros(nel,3);
                            
  fprintf('computing P1 flux jumps... ')
                                               
% Loop over Gaussian points
  for ig = 1:ngpt
    
      sigpt = oneg(ig);
      wt = onew(ig);
      
      % Compute flux jumps
      sigpt_ref = (1.0 + sigpt)/2.0;  % Map from [-1,1] to [0,1]
      [jmp] = p1fluxjmps_fnhandle(p1sol,eex,xy,evt,eboundt,tve,sigpt_ref,diff);

      % Loop over only the 3 edge bubble functions
      for j = 1:3
          % Original lines:
          [s,t] = reorder_s_x(sigpt,evt,j);   [psi_v,~,~] = vtqderiv(s,t,xl_v,yl_v);
          % LR modification:
          %[s,t] = reorder_s_leo(sigpt,j); %[psi_v,dpsidx_v,dpsidy_v] = tqderiv(s,t,xl_v,yl_v);
          edgeres(:,j) = edgeres(:,j) + (1/2) * wt * jmp(:,j) .* psi_v(:,j+3) .* (els(:,j)./2);
      end  
      % end edge bubble functions loop
    
  end
% end loop over Gaussian points
   
  fprintf('done\n')

end   % end function
