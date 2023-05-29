filename_in = '/home/ben/Documents/GitHub/PhD/third_year/paper/latex/EfficientAdaptiveSCforAdvDiff/figures/rf/stddev_tab4.dat';
filename_out = '/home/ben/Documents/GitHub/PhD/third_year/paper/latex/EfficientAdaptiveSCforAdvDiff/figures/rf/stddev_tab4_new.dat';

in = readmatrix(filename_in);

% Original mesh
x_in = in(:,1);
y_in = in(:,2);
u_in = in(:,3);

n=length(x_in);
nsqrt = sqrt(n);
% Reshape
x_grid = reshape(x_in,[nsqrt,nsqrt]).';
y_grid = reshape(y_in,[nsqrt,nsqrt]).';
u_grid = reshape(u_in,[nsqrt,nsqrt]).';

% Finer x y
x_coarse = x_grid(1,:).';
y_coarse = y_grid(:,1);

x_fine = interp1(1:nsqrt,x_coarse,1:0.5:nsqrt);
y_fine = interp1(1:nsqrt,y_coarse,1:0.5:nsqrt);

[X,Y] = meshgrid(x_fine,y_fine);


U = interp2(x_grid,y_grid,u_grid, X,Y);

Utr = U.'; Xtr = X.'; Ytr = Y.';

writematrix([Xtr(:), Ytr(:), Utr(:)], filename_out);
