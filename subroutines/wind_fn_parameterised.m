function w = wind_fn_parameterised(x,y,offset, factor)
xmod = factor*(x-offset(1)); ymod=factor*(y-offset(2));
zero_ind = xmod < -1 | xmod > 1 | ymod < -1 | ymod > 1;
wind_fn0 = @(x,y,nel) [2*y.*(1-x.*x), -2*x.*(1-y.*y)];
w = wind_fn0(xmod,ymod);
w(zero_ind,:) = 0;
end