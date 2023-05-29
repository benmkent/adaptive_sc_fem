function [u_interp] = scattered_interpolant_bk(connectivity, xy, U, xx_target, yy_target)

% Get values at vertices
for ivtx = 1:3
    xl_v(:,ivtx) = xy(connectivity(:,ivtx),1);
    yl_v(:,ivtx) = xy(connectivity(:,ivtx),2);
    ul_v(:,ivtx) = U(connectivity(:,ivtx));
end

% Form the triangulation
tri = triangulation(connectivity, xy);

% convert cartesian to barycentric
for ii = 1:length(xx_target)
    tri_id(ii) = pointLocation(tri,xx_target(ii),yy_target(ii));
    bary(ii,:) = cartesianToBarycentric(tri,tri_id(ii),[xx_target(ii),yy_target(ii)]);
    u_interp(ii) = sum(ul_v(tri_id(ii),:) .*bary(ii,:));
end