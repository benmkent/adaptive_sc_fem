function [DT,U,Udot,time,n ,nrej, dt, UDD] = tr_modified(A,M,f,uzero,dtzero, t0, tfinal,tol,nstar,info, udd_old)
%STABTR  stabilised TR integrator for n-dimensional system of ODEs
%   [DT,U,Udot,time] = stabtr(A,M,f,uzero,dtzero,tfinal,tol,10,0);
%   input
%          A, M      specified ODE system: M udot + A u = f without BCs
%          uzero     initial condition
%          dtzero    initial timestep
%          tfinal    final time
%          tol       local accuracy error tolerance
%          nstar     averaging frequency
%          info      output switch
%   output
%          DT        timestep history
%          U         solution history
%          Udot      solution time derivative history
%          time      discrete time evolution
%   X-IFISS function: DJS; 15 March 2018
% Copyright (c) 2018 D.J. Silvester, Huining Yang

%--------- setup
nmax=1e10;  % maximum number of timesteps
ub=uzero; nd=length(ub);
T=tfinal; dt0=dtzero;
%fprintf('StabTR iteration in %g dimensions ...',nd)
% preallocate array dimensions
DT = zeros(1,10); U = zeros(nd,10); time = zeros(1,10); Udot = zeros(nd,10);
UDD = zeros(nd,10);

%--------- initialization
udotb = M\(-A*ub+f(t0));                                % du/dt(0)
ke = sqrt((ub'*(M*ub))); acc = sqrt((udotb'*(M*udotb)));
itntable(nd,0,t0,0,ub,udotb,ke,acc)

%--------- first time step
v = (M+ 0.5*dt0*A)\(f(t0+dt0) + M*udotb - A*ub);               % first TR step
% v = (M+ 0.5*dt0*A)\(M*ub + 0.5*dt0*(f(t0) + f(t0+dt0) - A*ub));               % first TR step
u = ub + 0.5*dt0 *v;
udot = 2*(u-ub)/dt0 - udotb;	                       % du/dt(dt0)
udd = (udot-udotb)/dt0;		                           % second derivative
t=t0 + dt0;
r = norm(M*udot+A*u-(f(t)));
ke = sqrt((u'*(M*u))); acc = sqrt((udot'*(M*udot)));
itntable(nd,t,r,0,u,udot,ke,acc)

if ~isnan(udd_old)
%     u_tr = u;
%     u_ab = ub + 1.5 * dt0 * udot - 0.5 *udotb;
%     udiff = u_tr - u_ab;
    w =  udot + 0.5*dt0*udd_old;
    udiff =  0.5*v-w;
    d = (dt0*dt0/(3*(dt0+dt0)))*sqrt(abs(udiff'*M*udiff));
    dt = dt0; % (dt0*(tol/d)^(1/3));
else
    udd_old = nan(nd,1);
    dt = dt0;
end

n=2;                                                   % set time step index
DT(1:2) = [dt0, dt];
U(:,1) = ub; U(:,2) = u;
Udot(:,1) = udotb; Udot(:,2) = udot;
time(1:2) = [t0,t0+dt0];
UDD(:,1:2) = [udd_old, udd];

flag = 0; nrej=0; avflag = 0;  nav=nstar; tstar=tfinal;

%--------- loop until time limit is reached
while t <= T  & flag==0
    %if t+dt>T
    %    dt = T-t; flag = 1;
    %end                   % fix final time step
    if n==nmax, flag=1;
        %fprintf('\nToo slow -- step limit reached!'),
    end

    %   if condest(M+0.5*dt*A) < 1e-18
    %       fprintf(1,'dt = %.3e\n', dt);
    %       pause(0.1);
    %   end

    v = (M+0.5*dt*A)\(f(t+dt) + M*udot - A*u);                 % general TR step
    w = udot + 0.5*dt*udd;                               % AB2 step
    udiff =  0.5*v-w;
    upTR  = u + 0.5*dt*v;
    d = (dt*dt/(3*(dt+dt0)))*sqrt(abs(udiff'*M*udiff));

    % timestep control
    %   if d < ((1/0.7)^3)*tol*dt      % time step accepted
    if d < ((1/0.7)^3)*tol     % time step accepted
        if ((t>tstar & avflag==0) | ~rem(n,nav)) && ~flag  %%% smooth by averaging
            dt0 = 0.5*(dt+dt0);
            ub  = 0.5*(u+ub);
            udotb = 0.5*(udot+udotb);
            u = 0.5*(u + upTR);
            udot = 0.5*v;
            t = t + 0.5*dt;
            avflag=1;
            if nav == 1e4, nav = n; end
            if info==1, disp(['Averaging: n = ' int2str(n) ', t = ',num2str(t)]), end
        else                                   %%% take regular step
            dt0 = dt; t = t+dt0;
            ub = u; u = upTR;
            udotb = udot; udot = v - udot;
        end
        udd = (udot-udotb)/dt0;
        r = norm(M*udot+ A*u-f(t));                      % sanity check
        n=n+1;
        if n > length(time)		                    % allocate more memory
            DT   = [DT   zeros(1,100)];
            U    = [U zeros(nd,100)];   Udot = [Udot zeros(nd,100)];
            time = [time zeros(1,100)];
        end
        if n==3
            udiff_out = udiff;
        end
        DT(n) = dt; U(:,n) = u;  Udot(:,n) = udot; time(n) = t;UDD(:,n) = udd;
    else   % rejected step
        nrej = nrej + 1;
        if info==1, disp(['oops .. step ', int2str(n),' rejected']), end
    end

    % compute the next timestep
    dt = dt0; %(dt*(tol/d)^(1/3));
    ke = sqrt((u'*(M*u))); acc = sqrt((udot'*(M*udot)));
    itntable(nd,t,r,d,u,udot,ke,acc)
end
%--------- end of timestep loop


%--------- finishing touches
%fprintf('finished in %3i steps!\n',n)
DT = DT(1:n); U = U(:,1:n); Udot = Udot(:,1:n); time = time(1:n); UDD = UDD(:,1:n);
if nrej>0,
    %disp([int2str(nrej),' rejections in total: tol = ',num2str(tol)]),
end

% BMK Calculate the next value of w to feed forward
% w = udot + 0.5*dt*udd;                               % AB2 step
% feed forward udotb to calculate next w
return


function itntable(nd,t,r,d,u,udot,ke,acc)
% if nd==1,  % one-dimensional problem
%    if t==0; fprintf('\n t          u           udot       d\n'),end
%    fprintf(' %6.2e  %9.4e  %9.4e  %9.4e\n',t,u(1),udot(1),d)
% end
% if nd>1,  % nd-dimensional problem
% if t==0; fprintf('\n     t        ke         acc          d          res  \n'), end
% fprintf(' %6.2e  %9.4e  %9.4e   %9.4e  %9.4e\n',t,ke,acc,d,r)
% end
return
