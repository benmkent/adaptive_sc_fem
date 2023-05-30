function z_struct_new = propagate_sc_point(z_struct, tplusdt, params, fem)
if tplusdt > z_struct.t_z
    D = z_struct.D;
    N = z_struct.N;
    Q = z_struct.Q;
    f = z_struct.f;
    bc = z_struct.bc;
    bc_prime = z_struct.bc_prime;


    bound = z_struct.bound;
    notbound = z_struct.notbound;

    uzero = z_struct.u_z(:,end);
    dtzero = z_struct.dt_z(end);
    t0 = z_struct.t_z(end);

    tfinal = tplusdt;
    tol = params.letol;
    nstar = intmax;
    info = 0;

    %% Define system after Dirichlet bc
    Qx = Q; Qx(bound,:) = [];
    Qx(:,bound) = [];
    A = D+N;
    Ax = A; Ax(bound,:) = [];
    Ax(:,bound) = [];
    S = speye(size(Q));
    S(bound,:) = [];
    fx = @(t) S *f - A(notbound,bound)*bc(t) - Q(notbound,bound)*bc_prime(t);
    uzerox = uzero(notbound);


    if isempty(z_struct.w)
        w = zeros(size(Ax,1),0);
    else
        w = z_struct.w(:,end);
    end

    switch params.timestepping
        case 'stabtr'
        % Use modified stabtr algorithm
        %     [DT,U,Udot,time, n ,nrej, dt, UDD] = stabtr_modified(A,M,f,uzero(fem.notbound),dtzero,t0, tfinal,tol,nstar,info, w);
        [DT,U,Udot,time, n ,nrej, dt, UDD] = stabtr_modified(Ax,Qx,fx,uzerox,dtzero,t0, tfinal,tol,nstar,info, w);
        
        case 'fixed'
            tol = inf;
            [DT,U,Udot,time, n ,nrej, dt, UDD] = tr_modified(Ax,Qx,fx,uzerox,dtzero,t0, tfinal,tol,nstar,info, w);
    end

        if z_struct.tplusdt ~= time(1);
            U_with_bc(notbound,:) = [z_struct.u_z_tplusdt(notbound,:), U];
            time = [z_struct.tplusdt, time]; % Must also include previous alignement time to ensure can interpolate on whole interval.
            DT = [z_struct.dt_z_tplusdt, DT];
            UDD = [nan(size(UDD,1),1), UDD];
        else
            U_with_bc(notbound,:) = [U];
        end

    U_with_bc(bound,:) = z_struct.bc(time);

    % Interpolate to give approx at time tplusdt;
    u_z_tplusdt = interp1(time, [U_with_bc(:,:)].',tplusdt).';
    dt_z_tplusdt = interp1(time, [DT],tplusdt);

    % Construct updated data structure
    z_struct_new = z_struct;
    z_struct_new.tplusdt = tplusdt;
    z_struct_new.u_z = U_with_bc(:,:);
    z_struct_new.t_z = time(:);
    z_struct_new.u_z_tplusdt = u_z_tplusdt;
    z_struct_new.dt_z_tplusdt = dt_z_tplusdt;
    z_struct_new.dt_z = DT(:);
    z_struct_new.w = UDD(:,:);
    z_struct_new.n_steps = z_struct.n_steps+n;

    z_struct.downsampled = false;
else
    z_struct_new = z_struct;
    % We have stored the entire timestep history on this interval [t,t+dt) so can interpolate.
    % Interpolate to give approx at time tplusdt;
    %         if tplusdt < z_struct.t_z(1)
    % I think this is impossible as t (aka old tplusdt) is t_z(1)
    %             u_z_tplusdt = interp1([z_struct.tplusdt, z_struct.t_z], [z_struct.u_z_tplusdt, z_struct.u_z].',tplusdt).';
    %         else
    u_z_tplusdt = interp1([z_struct.t_z], [z_struct.u_z].',tplusdt).';
    dt_z_tplusdt = interp1([z_struct.t_z], [z_struct.dt_z],tplusdt);
    %         end
    if tplusdt < min([z_struct.t_z])
        error('Time history is not long enough');
    end
end
z_struct_new.u_z_tplusdt = u_z_tplusdt;
z_struct_new.tplusdt = tplusdt;
z_struct_new.dt_z_tplusdt = dt_z_tplusdt;
end