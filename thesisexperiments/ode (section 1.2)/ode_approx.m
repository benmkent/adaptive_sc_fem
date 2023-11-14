u_0 = 1;
diffcnst = 0.1;
a=1;

t = logspace(-3,2,1e3);
u = @(t,y) u_0 .* exp(-diffcnst * t)' .* exp(1i * t' *y);
exp_u = @(t) u_0 .* exp(-diffcnst * t(:)) .* sinc(t(:)/pi);
var_u = @(t) u_0^2 .* exp(-2* diffcnst * t(:)) .* (1- (sinc(t(:)/pi)).^2);
close all
figure(1); hold on;
plot(t,exp_u(t),'DisplayName','true');
figure(3); hold on;
plot(t,sqrt(var_u(t)),'DisplayName','true');

figure(7);
plot(t,sqrt(var_u(t)),'DisplayName','true')
set(gca,'XScale','log');

N = 1;
Sref = smolyak_grid_quick_preset(N,11-1);
SrRef = reduce_sparse_grid(Sref);

results =[0;t(:)];

%% Generate Colloc Pts
for w = 1:6;
N=1; % weights are already scaled for probability
S = smolyak_grid_quick_preset(N,w-1);
Sr = reduce_sparse_grid(S);
S2 = smolyak_grid_quick_preset(N,w-1+1);
Sr2 = reduce_sparse_grid(S2);

u_z = @(t) u(t,Sr.knots);

u_z_interp = interpolate_on_sparse_grid(S,Sr,u_z(t),Sr2.knots);
u_ref_interp = interpolate_on_sparse_grid(S,Sr,u_z(t),SrRef.knots);

exp_z = u_z(t) * Sr.weights';
var_z = (abs(u_z_interp)).^2 * Sr2.weights' - exp_z.^2;
std_z = sqrt((abs(u_z_interp)).^2 * Sr2.weights' - exp_z.^2);
e_est = sqrt((abs(u(t,Sr2.knots) - u_z_interp)).^2 * Sr2.weights');
e_ref = sqrt((abs(u(t,SrRef.knots) - u_ref_interp)).^2 * SrRef.weights');
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
% figure(8); hold on;
% plot(t, abs(sqrt(var_u(t)) - sqrt(var_z)'),'DisplayName',num2str(w));
% set(gca,'XScale','log','YScale','log');
% figure(5); hold on;
% plot(t,e_est,'DisplayName',num2str(w));
% set(gca,'XScale','log','YScale','log');
% 
% figure(6); hold on;
% plot(t,abs((e_ref - e_est)./e_ref),'DisplayName',num2str(w));
% set(gca,'XScale','log','YScale','log');

results = [results, [length(Sr.knots); e_ref(:)]];

end
results = [results, [length(SrRef.knots); 0*e_ref(:)]];
writematrix(results,'ode-approx-error.dat')
