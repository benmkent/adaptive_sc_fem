function [elerr] = flatten_spatial_estimator(eta_x_z_T,I_r)
%FLATTEN_SPATIAL_ESTIMATOR Combines spatial estimators across all
%collocation points
%
%   Inputs eta_x_z_T    Element errors for all collocation points
%          I_r         Reduced sparse grid
%   Outputs  elerr       Single element error representing all colloc points

% Combine element error estimators
elerr_z = zeros(length([eta_x_z_T{1}]),I_r.size);
for ii = 1:I_r.size
    elerr_z(:,ii) = [eta_x_z_T{ii}];
end

%% Use max error across all collocation points.
elerr = max(elerr_z,[],2);