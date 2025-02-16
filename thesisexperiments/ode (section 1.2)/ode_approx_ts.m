clear all;
N=1;
u_0 = 1;
diffcnst = 0.1;
a=1;

t = logspace(-3,3,200);
u = @(t,y) u_0 .* exp(-diffcnst * t)' .* exp(1i * t' *y);
exp_u = @(t) u_0 .* exp(-diffcnst * t) .* sinc(t/pi);
var_u = @(t) u_0^2 .* exp(-2* diffcnst * t) .* (1- (sinc(t/pi)).^2);
close all
% figure(1); hold on;
% plot(t,exp_u(t),'DisplayName','true');
% figure(3); hold on;
% plot(t,var_u(t),'DisplayName','true');

N=1;
Sref = smolyak_grid_quick_preset(N,10);
SrRef = reduce_sparse_grid(Sref);

%% Generate Colloc Pts
% Weights are already scaled, MIse not as defined in paper
w = 1;
S = smolyak_grid_quick_preset(N,w-1);
Sr = reduce_sparse_grid(S);
S2 = smolyak_grid_quick_preset(N,w+1);
Sr2 = reduce_sparse_grid(S2);


for ii = 1:Sr.n
    y = Sr.knots(ii);

    dt=1e-1; T=t(end);
    [u_tr, t_tr] = tr_bk(u_0,dt,T,y,diffcnst,a);

    A = @(y) -(-diffcnst + 1i * y);
    M=1;
    f = 0;
    uzero=u_0;
    dtzero=dt*1e-2;
    tfinal=T;
    tol=1e-2;
    nstar=inf;
    info=0;
    [DT,U,Udot,time]  = stabtr_clean(A(y),M,f,uzero,dtzero,tfinal,tol,nstar,info);

    error_tr = abs(u_tr - u(t_tr',y));
    error_stabtr = abs(U - u(time',y));

    figure(1); hold on;
    plot(t_tr, error_tr,'DisplayName',['Steps=' num2str(length(t_tr)-1) ', y=' num2str(y)]);
    set(gca,'XScale','log','YScale','log');

    figure(2); hold on;
    plot(time, error_stabtr,'DisplayName',['Steps=' num2str(length(time)-1) ', y=' num2str(y)]);
    set(gca,'XScale','log','YScale','log');

    figure(3); hold on;
    plot(t_tr,dt* ones(size(t_tr)),'DisplayName',['Steps=' num2str(length(t_tr)-1) ', y=' num2str(y)]);
    set(gca,'XScale','log','YScale','log');

    figure(4); hold on;
    plot(time,DT,'DisplayName',['Steps=' num2str(length(time)-1) ', y=' num2str(y)]);
    set(gca,'XScale','log','YScale','log');
end
%
%
% u_z_interp = interpolate_on_sparse_grid(S,Sr,u_z(t),Sr2.knots);
% u_ref_interp = interpolate_on_sparse_grid(S,Sr,u_z(t),SrRef.knots);
%
% exp_z = u_z(t) * Sr.weights';
% var_z = (abs(u_z_interp)).^2 * Sr2.weights' - exp_z.^2;
% e_est = sqrt((abs(u(t,Sr2.knots) - u_z_interp)).^2 * Sr2.weights');
% e_ref = sqrt((abs(u(t,SrRef.knots) - u_ref_interp)).^2 * SrRef.weights');
% figure(1);
% plot(t,exp_z,'DisplayName',num2str(w));
% set(gca,'XScale','log');
% figure(2); hold on;
% plot(t,abs(exp_z-exp_u(t)'),'DisplayName',num2str(w));
% set(gca,'XScale','log','YScale','log');
% figure(3);
% plot(t,var_z,'DisplayName',num2str(w))
% set(gca,'XScale','log');
% figure(4); hold on;
% plot(t, abs(var_u(t) - var_z'),'DisplayName',num2str(w));
% set(gca,'XScale','log','YScale','log');
% figure(5); hold on;
% plot(t,e_est,'DisplayName',num2str(w));
% set(gca,'XScale','log','YScale','log');
%
% figure(6); hold on;
% plot(t,abs((e_ref - e_est)./e_ref),'DisplayName',num2str(w));
% set(gca,'XScale','log','YScale','log');
%
% end
