function plot_wind_fn(wind_fn)

nplot = 1e1;
x_plot = linspace(-1,1,nplot);

y = rand(length(wind_fn),1);
y(1) = 1;

[X,Y] = meshgrid(x_plot);

for ii = 1:nplot
    for jj = 1:nplot
        w = [0,0];
        for kk = 1:length(wind_fn)
            w=w+ y(kk) * wind_fn{kk}(X(ii,jj),Y(ii,jj),1);
        end
        Z(ii,jj,1) = w(1);
        Z(ii,jj,2) = w(2);
    end
end

figure(5);
quiver(X,Y,Z(:,:,1),Z(:,:,2));
