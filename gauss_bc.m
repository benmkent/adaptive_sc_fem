function [ff] = gauss_bc(s,t,xl,yl,fem,time)
nel=length(xl(:,1));
zero_v = zeros(nel,1); xx=zero_v; yy=xx;
[xi,dxids,dxidt] = tshape(s,t);
for ivtx=1:3
    xx = xx + xi(ivtx) * xl(:,ivtx);
    yy = yy + xi(ivtx) * yl(:,ivtx);
end
ff=-specific_rhs_prime(xx,yy,nel,fem,time) - ;
ff = ff(:);
return

function f = specific_rhs(xx,yy,nel,fem,t)
BC(fem.notbound) = 0;
BC(fem.bound) = fem.bc_fn(t);
f = scattered_interpolant_bk(fem.T, fem.xy, BC, xx,yy);

function f = specific_rhs_prime(xx,yy,nel,fem,t)
BC(fem.notbound) = 0;
BC(fem.bound) = fem.bc_fn_prime(t);
f = scattered_interpolant_bk(fem.T, fem.xy, BC, xx,yy);