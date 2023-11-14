u_0 = 1;
diffcnst = 0.1;
a=1;

t = logspace(-3,2,200);
exp_u = u_0 .* exp(-diffcnst * t) .* sinc(t);
var_u = u_0^2 .* exp(-2* diffcnst * t) .* (1- (sinc(t)).^2);

figure(1)
plot(t,exp_u);
set(gca,'XScale','log');
figure(2)
plot(t,var_u);
set(gca,'XScale','log');
