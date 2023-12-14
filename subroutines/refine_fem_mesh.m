function [Z_I_star,fem_new] = refine_fem_mesh(Z_I_star,Mset,fem,params)
%REFINE_FEM_MESH Refines spatial meshes for all collocation points
% Inputs        Z_I_star    Collocation point data structure
%               Mset        Marked elements
%               fem         FEM data structure
%               params      Approximaton parameters
% Outputs       Z_I_star    Updated data structure
%               fem_new     Updated FEM structure

%% Set up
fem_new = fem;
% Construct detail space
[evtY,xyY,boundY,Ybasis] = p1grid_detail_space(fem.xy,fem.ev);
% Get marked elements and edges
[MMele,MMedge] = get_all_marked_elem(Ybasis,evtY,Mset,1);

% Refine mesh!
fprintf('Mesh refinement...');
meshRefTime = tic;
% Mesh refinement
[evt_new,xy_new,bound_new,interior_new,eboundt_new] = mesh_ref_bk(MMele,MMedge,fem.ev,fem.xy,fem.bound,evtY,xyY,boundY);
fprintf('done (%.5f sec)\n',toc(meshRefTime));

%% Redefine mass matrix and mesh properties
switch params.grid
    case'q1'
        [Q] = femq1_mass(xy,ev);

        [hx,hy,eex] = edgegen(xy,ev);
        tve=nan;
    case 'p1'
        [Q] = femp1_mass(xy_new,evt_new);

        [eex,tve,hx] = tedgegen(xy_new,evt_new);
        hy = 0*hx;

        T = evt_new;
        fem_new.T = T;
    otherwise
        error('Only P1 or Q1');
end

notbound = 1:size(xy_new,1);
notbound = setdiff(notbound, bound_new);

fem_new.xy = xy_new;
fem_new.ev = evt_new;
fem_new.bound = bound_new;
fem_new.notbound = notbound;
fem_new.hx = hx;
fem_new.hy = hy;
fem_new.eex = eex;
fem_new.ebound = eboundt_new;
fem_new.tve = tve;
fem_new.Q = Q;

%% Interpolate each approximation onto the new mesh
for ii = 1:length(Z_I_star)
    t(ii) = Z_I_star{ii}.t_z(1);
    U(:,ii) = Z_I_star{ii}.u_z(:,1);
    U_remesh(:,ii) = scattered_interpolant_bk(fem.T, fem.xy, U(:,ii), fem_new.xy(:,1),fem_new.xy(:,2));

    Z_I_star{ii}.u_z = U_remesh(:,ii);
    Z_I_star{ii}.u_z_tplusdt = U_remesh(:,ii);
    Z_I_star{ii}.t_z = t(ii);
    Z_I_star{ii}.tplusdt = t(ii);
    Z_I_star{ii}.w = nan(1);
    Z_I_star{ii}.bound = bound_new;
    Z_I_star{ii}.notbound = notbound;
end

%% Plot the new spatial mesh
if 1==1 % 
% plotting the generated grid
    nvtx=length(fem_new.xy(:,1));
    rnel=length(fem_new.ev(:,1));
	adj=sparse(nvtx,nvtx);
    for i=1:rnel
	adj(fem_new.ev(i,1),fem_new.ev(i,2)) =1;
	adj(fem_new.ev(i,2),fem_new.ev(i,3)) =1;
	adj(fem_new.ev(i,3),fem_new.ev(i,1)) =1;
    end
    figure(10)
    gplot(adj,fem_new.xy,'b')
    axis('square')
    hold on
    plot(fem_new.xy(:,1),fem_new.xy(:,2),'ro')
    rxybd=fem_new.xy(fem_new.bound,:);
    plot(rxybd(:,1),rxybd(:,2),'ko')
    hold off
    title('P1 finite element subdivision');
    drawnow
end