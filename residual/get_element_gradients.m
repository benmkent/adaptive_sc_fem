function [dUx, dUy] = get_element_gradients(p1sol, fem)
%GET_ELEMENT_GRADIENTS Computes the spatial gradients for a FE
%approximation
%
% Inputs    p1sol   approximation on p1 mesh
%           fem     finite element data structure
% Outputs   dUx     gradient in spatial dimension 1 on each element
%           dUy     gradient in spatial dimension 2 on each element

%% Set up
xy = fem.xy;
evt = fem.ev;
x = xy(:,1);
y = xy(:,2);
nel = length(evt(:,1));

%% Recover local coordinates and local solution
for ivtx = 1:3
    xl_v(:,ivtx) = x(evt(:,ivtx));
    yl_v(:,ivtx) = y(evt(:,ivtx));
    sl_v(:,ivtx) = p1sol(evt(:,ivtx));
end

%% For each element get normal vector of sln plane.
% P1 approximation implies solution is a plane and derivatives
% are constant on the element.
v1(:,1:3) = [xl_v(:,2) - xl_v(:,1),...
    yl_v(:,2) - yl_v(:,1),...
    sl_v(:,2) - sl_v(:,1)];
v2(:,1:3) = [xl_v(:,3) - xl_v(:,1),...
    yl_v(:,3) - yl_v(:,1),...
    sl_v(:,3) - sl_v(:,1)];

normvectors = cross(v1,v2,2);

%% Compute gradients
% Each plane is of the form nx x + ny y + nu u = 0
% therefore u = -(nx/nu) x -(ny/nu) y
% and partial derivatives are simply
dUx = -normvectors(:,1)./normvectors(:,3);
dUy = -normvectors(:,2)./normvectors(:,3);