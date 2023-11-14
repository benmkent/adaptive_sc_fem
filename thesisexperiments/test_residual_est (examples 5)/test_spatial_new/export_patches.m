
level = 'l4';
% timeInd = 6; timeStr='t01';
% timeInd = 16; timeStr='t1';
% timeInd = 29; timeStr='t10';
timeInd = 42; timeStr='t100';

xy = fem.xy;
evt = fem.ev;

sol = data_table{timeInd,'u_z'}{1};
err = data_table{timeInd,'pi_x_z_T'}{1};
name = [level timeStr];

%% Recover local coordinates of elements
nel = size(evt,1);
xl_v = zeros(nel,3);
yl_v = zeros(nel,3);
for ivtx = 1:3
    xl_v(:,ivtx) = xy(evt(:,ivtx),1); % x-coordinates of all vertices
    yl_v(:,ivtx) = xy(evt(:,ivtx),2); % y-coordinates of all vertices
    sol_v(:,ivtx) = sol(evt(:,ivtx));
end

%% Write patch file
patches_sol = zeros(0,3);
for ii= 1:nel
    patches_sol = [patches_sol; xl_v(ii,:).',yl_v(ii,:).',sol_v(ii,:).'];
end

patches_err = zeros(0,3);
for ii= 1:nel
    patches_err = [patches_err; xl_v(ii,:).',yl_v(ii,:).',err(ii)*[1;1;1]];
end

writematrix(patches_sol,[name 'patches_sol.dat'],'Delimiter','space');
writecell([{'x','y','c'}; num2cell(patches_err)],[name 'patches_err.dat'],'Delimiter','space');

%writematrix(evt-1,'evt.dat','Delimiter','space')
%writematrix([xy, sol],'sol.dat','Delimiter','space')