function [u_xy, x,y] = vec2xy(u_vec, fem)
n_x = sqrt(size(fem.xy,1));
xy=fem.xy;
x = reshape(xy(:,1), [n_x,n_x]);
y = reshape(xy(:,2), [n_x,n_x]);

u_xy = reshape(u_vec, [n_x,n_x]);
end