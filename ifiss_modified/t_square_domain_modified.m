function [ mv, xy, bound, mbound, grid_type, outbc, x, y] = square_domain(square_type, nc, grid_type)
%SQUARE_DOMAIN Modified TIFISS function   
% square domain Q2 grid generator
%   Modified BMK 2022: no longers saves out to function or requests user
%   inputs
%
%   square_domain(square_type, grid_reg);
%   input
%          square_type   1 for  [0,1]x[0,1]
%                        2 for [-1,1]x[-1,1]               
%          grid_reg      1 for uniform grid
%                        0 for option to stretch grid              
% 
% grid defining data is saved to the file: square_grid.mat
%    PIFISS function: DJS; 5 February 2007. 
% Copyright (c) 2007 C.E. Powell, D.J. Silvester
if square_type==1,
   fprintf('\n\nGrid generation for unit square  domain.\n')
   elseif square_type==2,
   fprintf('\n\nGrid generation for reference square domain.\n')
   else
   error('illegal input parameter: square_type')
end
% nc=default('grid parameter: 3 for underlying 8x8 grid (default is 16x16)',4);
if nc<2, error('illegal parameter choice, try again.'), end
% if grid_reg==0
%    grid_type=default('uniform/stretched grid (1/2) (default uniform)',1);
% else, grid_type=1; end
n=2^nc; np=n/2; nq=n/4;
%
%% compute (x,y) coordinates of vertices
% y-direction
if grid_type==2
   hmax=nc/(2^(nc+1));
   x1=-1;x2=-2*hmax;x3=2*hmax;x4=1;nx1=2^(nc-1)-1;nx2=2;nx3=2^(nc-1)-1;
   y1=-1;y2=-2*hmax;y3=2*hmax;y4=1;ny1=2^(nc-1)-1;ny2=2;ny3=2^(nc-1)-1;
   y=subint(y1,y2,y3,y4,ny1,ny2,ny3);
   stretch=(y(3)-y(2))/(y(2)-y(1));
   left=-1;
   x=y;
else
   yy=[1/np:1/np:1];
   ypos=[0,yy];
   yneg=-yy(length(yy):-1:1);
   y=[yneg,ypos]'; 
   left=-1;
   if square_type==1
      y=[0:1/(2*np):1]'; left=0;
   end
   x=y; 
   end
%
%% compute biquadratic element coordinates
nvtx=(n+1)*(n+1);
[X,Y]=meshgrid(x,y);
xx=reshape(X',nvtx,1);
yy=reshape(Y',nvtx,1);
xy=[xx(:),yy(:)];
%
kx = 1;
ky = 1;
mel=0;
for j=1:np
   for i=1:np
      mref=(n+1)*(ky-1)+kx;
      mel=mel+1;
      nvv(1) = mref;
      nvv(2) = mref+2;
      nvv(3) = mref+2*n+4;
      nvv(4) = mref+2*n+2;
      nvv(5) = mref+1;
      nvv(6) = mref+n+3; 
      nvv(7) = mref+2*n+3; 
      nvv(8)=  mref+n+1;
      nvv(9)=  mref+n+2; 
      mv(mel,1:9)=nvv(1:9);
      kx = kx + 2;
   end
   ky = ky + 2; 
   kx = 1;
end
%
%% compute boundary vertices and edges
% four boundary edges 
k1=find( xy(:,2)==left );
e1=[]; for k=1:mel, if any(mv(k,5)==k1), e1=[e1,k]; end, end
ef1=ones(size(e1));
%
k2=find( xy(:,1)==1  & xy(:,2)<=1   & xy(:,2) >left);
e2=[]; for k=1:mel, if any(mv(k,6)==k2), e2=[e2,k]; end, end
ef2=2*ones(size(e2));
%
k3=find( xy(:,2)==1  & xy(:,1)<1   & xy(:,1) >left);
e3=[]; for k=1:mel, if any(mv(k,7)==k3), e3=[e3,k]; end, end
ef3=3*ones(size(e3));
%
k4=find( xy(:,1)==left & xy(:,2)<=1   & xy(:,2) >left );
e4=[]; for k=1:mel, if any(mv(k,8)==k4), e4=[e4,k]; end, end
ef4=4*ones(size(e4));
%

bound=sort([k1;k2;k3;k4]);
mbound=[e1',ef1';e2',ef2';e3',ef3';e4',ef4'];
%
%
outbc=1;

% save square_grid.mat mv xy bound mbound grid_type outbc x y
 
return