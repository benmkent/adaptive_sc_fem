function plot_data(pmethod,dom_type,sol,errelem,evt,xy)
%PLOT_DATA plots solution as well as the error estimate
%
% plot_data(pmethod,dom_type,sol,errelem,evt,xy)
%
% input: 
%         pmethod      approximation method
%        dom_type      domain type
%             sol      nodal FE solution vector
%         errelem      element error indicators
%             evt      element mapping matrix
%              xy      vertex coordinate vector  
%
% NOTE that the solution and element-indicators are interpolated on a square 
% grid [X,Y] in order to plot isolines (contour plot); Matlab does not provide 
% a countour function for mesh-based functions.
%
%   TIFISS function: LR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, L. Rocchi

  global slope1crack 
  global slope2crack

% Recover correct element mapping matrix if P2 approximation is used  
  if isequal(pmethod,2)
     evt = evt(:,1:3);
  end
  
  nel  = size(evt,1);   % Number of elements

% Refine grid and get the cartesian product mesh
  npoints = 100; 
  x = linspace(min(xy(:,1)),max(xy(:,1)),npoints); 
  y = linspace(min(xy(:,2)),max(xy(:,2)),npoints); 
  [X1,Y1] = meshgrid(x,y);

% Gridded solution for contour plot 
  solcont = griddata(xy(:,1),xy(:,2),sol,X1,Y1); 
  
% Recover local coordinates of elements
  xl_v = zeros(nel,3); 
  yl_v = zeros(nel,3); 
  for ivtx = 1:3
      xl_v(:,ivtx) = xy(evt(:,ivtx),1); % x-coordinates of all vertices
      yl_v(:,ivtx) = xy(evt(:,ivtx),2); % y-coordinates of all vertices
  end
% Recover the element's centroid coordinates
  xyc(:,1) = (1/3) * sum(xl_v,2);
  xyc(:,2) = (1/3) * sum(yl_v,2);

% Refined grid and get cartesian product mesh
  x = 0.5 * ( x(1:end-1) + x(2:end) );
  y = 0.5 * ( y(1:end-1) + y(2:end) );
  [X2,Y2] = meshgrid(x,y); 
  
% Error estimator
  err = griddata(xyc(:,1),xyc(:,2),errelem,X2,Y2);
% Fix to zero (eventual very small) negative values (due to griddata interpolation)
  err(err<0.0) = 0.0;  
  
% Eliminate data outside the spatial domain  
  if dom_type==2 
      % L-shaped domain
      solcont(X1<0 & Y1<0) = nan;
      err(X2<0 & Y2<0)     = nan;
  elseif dom_type==3
      % Crack domain
      solcont((Y1>slope1crack*X1 & X1<0) & (Y1<0)) = nan;
      solcont((Y1<slope2crack*X1 & X1<0) & (Y1>0)) = nan;
      solcont((Y1>slope1crack*X1 & X1<0) & (Y1<0)) = nan;
      solcont((Y1<slope2crack*X1 & X1<0) & (Y1>0)) = nan;
      err((Y2>slope1crack*X2 & X2<0) & (Y2<0))     = nan;    
      err((Y2<slope2crack*X2 & X2<0) & (Y2>0))     = nan;
      err((Y2>slope1crack*X2 & X2<0) & (Y2<0))     = nan;    
      err((Y2<slope2crack*X2 & X2<0) & (Y2>0))     = nan;
  end  
  
% -----------------------------------------------------------------------------  
% Plot solution and spatial estimate
% -----------------------------------------------------------------------------   
  figure;
  
  subplot(221)
  contour(X1,Y1,solcont,20);
  axis square;  axis off; 
  title(['Solution with P',num2str(pmethod),' approximations']);  %title(titleSolution);
  outlinedomain(dom_type,xy);
  
  subplot(222)
  trimesh(evt,xy(:,1),xy(:,2),sol);
  axis square; view(330,30);
  
  subplot(223)
  contour(X2,Y2,err,20);
  axis square;  axis off; 
  title('Estimated error');
  outlinedomain(dom_type,xy);
 
  subplot(224)
  mesh(X2,Y2,err);           
  axis square; view(330,30);
  
  set(findall(gcf,'-property','Fontsize'),'Fontsize',13);  
end %end function
  

% -----------------------------------------------------------------------------
% Child function
% -----------------------------------------------------------------------------  
function outlinedomain(dom_type,xy)
% Calling the function outlining the domain in contour plot 
  if dom_type == 1
      if min(xy(:,1)) == 0
          unitsquare;
      else%min(xy(:,1)) == -1
          squarex;
      end
  elseif dom_type==2
      ellx;
  else%dom_type==3
      largecrack;
  end
end % end child function