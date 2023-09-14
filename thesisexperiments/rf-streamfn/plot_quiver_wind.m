load wind.mat
load approx.mat

fem.xy = xy;
params.grid = 'q1';

%%
% rng(0);
ii_z = randi(1371,1);
% z = knots5(:,ii_z);
z = knots100(:,ii_z);

figure(2);
plot(z);
% z(64)=1;
% z(1)=0.1;
wxz = tensorprod(wx,z,3,1);
wyz = tensorprod(wy,z,3,1);

wx0 = 2*(1-xp.^2).*yp;
wy0 = -2*(1-yp.^2).*xp;

figure(1); hold on;
ax = gca();
cla(ax);
params.plottype = 'contourf';
params.plotview=[0,90];
% plot_surf(ax,u5(:,ii_z),fem,params); view(0,90); axis square;
plot_surf(ax,u100(:,ii_z),fem,params); view(0,90); axis square;

[X,Y] = meshgrid(linspace(-1,1,16));
wxplot = interp2(xp.',yp.',(wx0+wxz).',X,Y);
wyplot = interp2(xp.',yp.',(wy0+wyz).',X,Y);
quiver(X,Y,wxplot,wyplot,'Color','white');
axis square
xlim([-1,1]);ylim([-1,1]);
view([0,90])