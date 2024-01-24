function fem = initialise_fem_mesh(params)
%INITIALISE_FEM_MESH Constructs the mesh objects using IFISS/TIFISS
% Input
%   params      Structure with user input mesh options
% Output
%   fem         Structure containing mesh properties.

%% Generate grid and mass matrix
% Switch for Q1 or P1
switch params.grid
    case 'q1'
        % See IFISS functions
        [mv, xy, bound, mbound, ~, ~, ~, ~] = ...
            rsquare_domain_modified(params.grid_param, params.grid_type);
        [ev,ebound] = q1grid(xy,mv,bound,mbound);
        [Q] = femq1_mass(xy,ev);
        [hx,hy,eex] = edgegen(xy,ev);
        tve=nan;
    case 'p1'
        % See TIFISS functions
        % Construct P1 mesh from initial level 2 mesh and subsequent
        % refinements.
        ini_grid_param = 2;
        [ mv, xy, bound, mbound, ~, ~, ~, ~] = t_square_domain_modified(2, ini_grid_param, params.grid_type);
        [ev,ebound] = p1grid(xy,mv,bound,mbound);
        for i=1:(params.grid_param-1)
            [xy,ev,bound,ebound] = p1_refinement(xy,ev,bound,ebound);
        end
        % Validate mesh
        [Q] = femp1_mass(xy,ev);
        [eex,tve,hx] = tedgegen(xy,ev);
        hy = 0*hx;

        T = ev;
        fem.T = T;
    otherwise
        error('Only P1 or Q1');
end

% Define indices that are not in boundary
notbound = 1:size(xy,1);
notbound = setdiff(notbound, bound);

%% Construct FEM structure holding mesh properties
fem.xy = xy;
fem.ev = ev;
fem.bound = bound;
fem.notbound = notbound;
fem.hx = hx;
fem.hy = hy;
fem.eex = eex;
fem.ebound = ebound;
fem.tve = tve;
fem.Q = Q;

fem.poincare = 2;
fem.amax=0.1;
fem.amin=0.1;