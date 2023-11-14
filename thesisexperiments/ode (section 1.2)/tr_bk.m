function [u_tr, t_tr] = tr_bk(u_0,dt,T,y,diffcnst,a)

theta = 0.5;
% dt = 1e-3;
u_tr(1) = u_0;
t_tr(1) = 0;
ii = 2;
while t_tr(end) < T
   u_tr(ii) = (1- dt*(1-theta)*(-diffcnst + 1i * y))^-1* (1 + dt* theta *(-diffcnst + 1i * y)) * u_tr(ii-1);
   t_tr(ii) = t_tr(ii-1) + dt;
   ii=ii+1;
end