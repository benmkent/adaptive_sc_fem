function z_struct = downsample_to_reference_times(z_struct, reference_times)
    if z_struct.downsampled == true
        %No processing required
    else
    t0 = z_struct.t_z(1);
    t1 = z_struct.t_z(end);
    
    t_in_t0_t1 = reference_times(reference_times > t0 & reference_times < t1);

    % We want to save the first approximation after the previous sync time,
    % i.e. u_z(2), and its corresponding data for propagting fwd in time.
    % This will allow the reconstruction of intermidiate times (before
    % t(end), by resolving the timestepping problem, but without having to
    % save all of the intermeditate timesteps which causes a memory issue.
    % This is clearly an additional computational cost.
    if length(z_struct.t_z) > 2
        t_z = [t0, z_struct.t_z(2), t_in_t0_t1, t1];
        [t_z, sort_indices] = sort(t_z);
         % Identify where the first timestep has been placed in the new
         % vector t_z. We need this as it is the only thing we can validly
         % interpolate from t0 to, and will provide data to propagate
         % forwards in time if we need some of the lost data.
        ind_tz2 = sort_indices(2);
    else
        t_z = [t0, t_in_t0_t1, t1];
        % If there were no other timesteps (i.e t_z is of length 2), then
        % only the final step is valid anyway so can interpolate linearly
        % between t0 and t1.
        ind_tz2 = length(t_z);
    end
    u_z = interp1(z_struct.t_z, z_struct.u_z',t_z)';
    dt_z = nan(size(t_z));
    dt_z(1) = z_struct.dt_z(1);
    dt_z(end) = z_struct.dt_z(end);
    w_z = nan(size(z_struct.w,1),length(t_z));
    w_z(:,1) = z_struct.w(:,1);
    w_z(:,end) = z_struct.w(:,end);
    
    % Save out the appropriate data for the index ind_tz2;
    dt_z(ind_tz2) = z_struct.dt_z(2);
    w_z(:,ind_tz2) = z_struct.w(:,2);

    z_struct.t_z = t_z;
    z_struct.u_z = u_z;
    z_struct.dt_z = dt_z;
    z_struct.w = w_z;

    z_struct.downsampled = true;
    end
end
