problem = define_problem('doubleglazing');
params = define_params('smolyak-ref');

params.grid_type = 1;
params.grid_param = 5;

problem.viscosity = 1e1;

fem = initialise_fem_matrices_new(problem, params);
y = zeros(4,1);
h = fem.hx(1);
Q = fem.mass;
D = fem.diff(y);
W = fem.conv(y);
H = fem.stiffness;

n = size(Q,1);

lambdaQ = eigs(Q,n);
lambdaD = eigs(D,n);
lambdaW = eigs(W,n);
lambdaH = eigs(H,n);
lambdaOp = eigs( D + W, -Q, 100,'smallestabs');

% figure(1); hold on;
% scatter(h, real(lambdaQ));
% set(gca,'XScale','log')
% set(gca,'YScale','log')
% 
% figure(2); hold on;
% scatter(h, lambdaD);
% set(gca,'XScale','log')
% set(gca,'YScale','log')
% 
% figure(3); hold on;
% scatter(h,imag(lambdaW));
% set(gca,'XScale','log')
% figure(4); hold on;
% scatter(-real(lambdaOp),imag(lambdaOp));
% set(gca,'XScale','log')
% scatter(real(lambdaQ),imag(lambdaQ));

writematrix([h*ones(n,1), lambdaQ],['lambdaQgridType' num2str(params.grid_type) 'gridParam' num2str(params.grid_param) '.dat']);
writematrix([h*ones(n,1),lambdaD],['lambdaDgridType' num2str(params.grid_type) 'gridParam' num2str(params.grid_param) '.dat']);
writematrix([h*ones(n,1),imag(lambdaW)],['lambdaWgridType' num2str(params.grid_type) 'gridParam' num2str(params.grid_param) '.dat']);
writematrix([h*ones(length(lambdaOp),1), -real(lambdaOp), imag(lambdaOp)],['lambdaOpgridType' num2str(params.grid_type) 'gridParam' num2str(params.grid_param) '.dat']);
