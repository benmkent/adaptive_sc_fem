function w = recirculating_parameterised(x,y,offset, factor)
%RECIRCULATING_PARAMETERISED Recirculating wind field with geometry changes
% Inputs
%   x       vector of x to evaluate
%   y       vector of y to evaluate
%   offset  pair(a,b) for centre of wind field
%   factor  scaling factor in x and y
% Outputs
%   w       wind field at (x,y) 

%% Define modified coordinates
xmod = factor*(x-offset(1)); ymod=factor*(y-offset(2));
zero_ind = xmod < -1 | xmod > 1 | ymod < -1 | ymod > 1;

% Evaluate recirculating wind field with modified (x,y) coordinates
wind_fn0 = @(x,y,nel) [2*y.*(1-x.*x), -2*x.*(1-y.*y)];
w = wind_fn0(xmod,ymod);
w(zero_ind,:) = 0;
end