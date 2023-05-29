function plot_datatable(data_table, reference,fem,params)
close all;
h(1) = figure(1);
%% Plot Exp StdDev
% Get xy grid
xy = fem.xy;
% Set 4 plot times
t_array = [0.1,1,10,100];
for ii = 1:4
%      Find index and select data
    ind_t = find([reference.data_table.t] >= t_array(ii) - eps,1);
    t = reference.data_table{ind_t,'t'};
    u = reference.data_table{ind_t,'u'};
    u = u{1};
    I = reference.data_table{ind_t,'I'};
    I = I{1};
%     Find index and select data
%     ind_t = find([reference.reference_table.t] >= t_array(ii) - eps,1);
%     t = reference.reference_table{ind_t,'t'};
%     u = reference.reference_table{ind_t,'u'};
%     u = u{1};
%     I = reference.reference_table{ind_t,'I'};
%     I = I{1};
    I_r = reduce_sparse_grid(I);
    C = get_mi_set(I);
    [~,C_star] = check_set_admissibility([C; sparse_grid_reduced_margin(C)]);
    C_star = unique(C_star,'rows');
    knot_fn = params.knot_fn;
    lev2knots = params.lev2knots;
    I_star = smolyak_grid_multiidx_set(C_star,knot_fn,lev2knots);
    I_star_r = reduce_sparse_grid(I_star);
    u_interpolated = interpolate_on_sparse_grid(I,I_r,u, I_star_r.knots);
%     
%     Calc exp and var
    exp_u = sum(u .* I_r.weights,2);
    exp_u2 = sum(u_interpolated.^2 .* I_star_r.weights,2);
    var_u = exp_u2 - exp_u.^2;
    var_u(var_u < 0) = 0;
    nx = sqrt(length(exp_u));
    [exp_u, x,y] = vec2xy(exp_u, fem);
    [var_u, x,y] = vec2xy(var_u, fem);
% 
%     Plot and write data
%     h(end+1)=figure;
%     contourf(x,y,exp_u);
%     title(num2str(t));
%     colorbar();
%     caxis([0,1]);
%     matlab2tikz('filename',['exp' num2str(ii) '.tex']);
    writematrix([x(:),y(:), exp_u(:)],['exp_tab' num2str(ii) '.dat'],"Delimiter"," ")
    writematrix([x(:),y(:), sqrt(var_u(:))],['stddev_tab' num2str(ii) '.dat'],'Delimiter','space')
% 
%     h(end+1)=figure;
%     contourf(x,y,sqrt(var_u));
%     title(num2str(t));
%     colorbar();
%     matlab2tikz('filename',['std' num2str(ii) '.tex']);
end

%% Plot data_table
% Select only time indices in which there is no refinement.
inds_valid = cellfun(@isempty,data_table{:,'pi_I_alpha'});

%% Plot error components
h(end+1) = figure();
plot(data_table{:,'t'},data_table{:,'E'},'LineStyle',':','DisplayName','Threshold')
hold on;
plot(data_table{:,'t'},data_table{:,'pi'},'DisplayName','Est')
plot(data_table{:,'t'},data_table{:,'pi_I'},'DisplayName','Interp')
plot(data_table{:,'t'},data_table{:,'pi_I_delta'},'DisplayName','Corr')
plot(data_table{:,'t'},data_table{:,'pi_delta'},'DisplayName','TS')
plot(data_table{:,'t'},data_table{:,'pi_delta'}+data_table{:,'pi_I_delta'},'DisplayName','Corr+TS')
plot(reference.data_table{:,'t'},reference.data_table{:,'error'},'DisplayName','True');
set(gca,'XScale','log','YScale','log');
legend('show');
hold off;
%matlab2tikz('filename',['components.tex']);
% writematrix([data_table{:,{'t','pi','pi_I','pi_I_delta','pi_delta'}}],'components-errorest.dat','Delimiter','space');
% writematrix([reference.data_table{:,{'t','error'}}],'error.dat','Delimiter','space');

%% Identify largest grid in use (final time of reference)
if isempty(reference.reference_table)
    I_super = data_table{end,'I_star'};
else
    I_super = reference.reference_table{end,'I'};
end
I_super = I_super{1};
I_super_r = reduce_sparse_grid(I_super);
ge_est_in_super = nan(size(data_table,1),I_super_r.size);
n_steps_in_super = zeros(size(data_table,1),I_super_r.size);
I_star_in_super = nan(size(data_table,1),I_super_r.size);
I_in_super = nan(size(data_table,1),I_super_r.size);
for ii = 1:size(data_table,1)
    I_star = data_table{ii,'I_star'};
    I_star = I_star{1};
    I_star_r = reduce_sparse_grid(I_star);
    if I_super_r.size == I_star_r.size
        inds_in_super = 1:I_super_r.size;
    else
        [~,inds_in_super,~,~] = compare_sparse_grids(I_super,I_super_r,I_star,I_star_r);
    end
    cell_ge_ii = data_table{ii,'ge_estimate'};
    ge_est_in_super(ii,inds_in_super) = cell_ge_ii{1};
    cell_nsteps_ii = data_table{ii,'n_steps'};
    n_steps_in_super(ii,inds_in_super) = cell_nsteps_ii{1};
    I_star_in_super(ii,inds_in_super) = 1;
    I = data_table{ii,'I'};
    I = I{1};
    I_r = reduce_sparse_grid(I);
    ii_zero_data(ii) = find(sum(abs(I_r.knots),1) < 1e-6,1,'first');

    [~,inds_in_super,~,~] = compare_sparse_grids(I_super,I_super_r,I,I_r);
    I_in_super(ii,inds_in_super) = 1;
    ge_zero(ii) = ge_est_in_super(ii,inds_in_super(ii_zero_data(ii)));
    n_steps_zero(ii) = n_steps_in_super(ii,inds_in_super(ii_zero_data(ii)));
%     J = data_table{ii,'J'};
%     J = J{1};
%     if ii==1
%         MI{1} = [ones(1,size(I_star_r.knots,1));J];
%     else
%         MI{ii} = [MI{ii-1}; J];
%     end
    MI = get_mi_set(I);
    MI_max(ii,:) = max(MI,[],1);
end
ge_ref_in_super = nan(size(reference.data_table,1),I_super_r.size);
Iref_in_super = nan(size(reference.data_table,1),I_super_r.size);
for ii = 1:size(reference.data_table,1)
    I_star = reference.data_table{ii,'I'};
    I_star = I_star{1};
    I_star_r = reduce_sparse_grid(I_star);
    ii_zero(ii) = find(sum(abs(I_star_r.knots),1) < 1e-6,1,'first');
    [~,inds_in_super,~,~] = compare_sparse_grids(I_super,I_super_r,I_star,I_star_r);
    cell_ge_ii = reference.data_table{ii,'ts_error'};
    ge_ref_in_super(ii,inds_in_super) = cell_ge_ii{1};
    Iref_in_super(ii,inds_in_super) = 1;
    tserror_zero(ii) = ge_ref_in_super(ii,inds_in_super(ii_zero(ii)));
end


%% Plot timestepping error estimates
h(end+1) = figure(); cla(); hold on; ax5 = gca();
% [e_mean,e_std, e_min, e_max] = plot_error_stats(ax5, reference.data_table{:,'t'}, reference.data_table{:,'ts_error'},[1,0,0]);
% [est_mean,est_std, est_min, est_max] = plot_error_stats(ax5, data_table{inds_valid,'t'}, data_table{inds_valid,'ge_estimate'},[0,1,0]);
plot(ax5, data_table{inds_valid,'t'}, mean(ge_est_in_super(inds_valid,:).*I_in_super(inds_valid,:),2,'omitnan'));
plot(ax5, data_table{inds_valid,'t'}, max(ge_est_in_super(inds_valid,:).*I_in_super(inds_valid,:),[],2,'omitnan'));
plot(ax5, reference.data_table{:,'t'}, mean(ge_ref_in_super.*Iref_in_super,2,'omitnan'));
plot(ax5, reference.data_table{:,'t'}, max(ge_ref_in_super.*Iref_in_super,[],2,'omitnan'));
plot(ax5, reference.data_table{:,'t'}, tserror_zero);
set(gca,'XScale','log','YScale','log');
%matlab2tikz('filename',['ge_error.tex']);

%% Plot error v estimators
est_mean = mean(ge_est_in_super.*I_in_super,2,'omitnan');
h(end+1) = figure(); cla(); hold on; ax = gca();
plot(data_table{:,'t'},data_table{:,'pi'});
plot(reference.data_table{:,'t'},reference.data_table{:,'error'},'Marker','x');
plot(data_table{:,'t'},data_table{:,'pi_I'} + est_mean);
legend('show'); set(gca,'XScale','log','YScale','log');
%matlab2tikz('filename',['error-est.tex']);

h(end+1) = figure(); cla(); hold on;
t_union = unique([reference.data_table{:,'t'};]); %data_table{:,'t'}]);

error_union = interp1(reference.data_table{:,'t'}, reference.data_table{:,'error'}, t_union,'linear');
est_union = interp1(data_table{inds_valid,'t'}, data_table{inds_valid,'pi'}, t_union,'linear');
plot(t_union, est_union./error_union,'x');
set(gca,'XScale','log','YScale','log');
est_union2 = interp1(data_table{inds_valid,'t'}, data_table{inds_valid,'pi_I'} + est_mean(inds_valid), t_union,'linear');
plot(t_union, est_union2./error_union,'x');
set(gca,'XScale','log','YScale','log');
%matlab2tikz('filename',['efficiency.tex']);

for ii = 1:size(I_in_super,2)
    error_ge_est_union(:,ii) = interp1(data_table{inds_valid,'t'}, ge_est_in_super(inds_valid,ii), t_union,'linear');
    efficiency_ge(:,ii) = error_ge_est_union(:,ii) ./ ge_ref_in_super(:,ii);
end
h(end+1) = figure(); axGEEfficiency = gca(); hold on;
plot(t_union, mean(efficiency_ge,2,'omitnan'));
plot(t_union, min(efficiency_ge,[],2,'omitnan'));
plot(t_union, max(efficiency_ge,[],2,'omitnan'));
set(axGEEfficiency,'XScale','log','YScale','log');

h(end+1) = figure(); cla();
plot(data_table{:,'t'},data_table{:,'delta_t'});
xlabel('t');ylabel('\Delta t^{alg}');
set(gca,'XScale','log','YScale','log');
%matlab2tikz('filename',['alg-timesteps.tex']);

h(end+1) = figure(); hold on; ax8 = gca(); cla();
[dt_mean,dt_std, dt_min, dt_max] = plot_error_stats(ax8, data_table{:,'t'}, data_table{:,'dt_z'},[1,0,0]);
xlabel('t'); ylabel('Timesteps');
set(ax8,'XScale','log','YScale','log');
%matlab2tikz('filename',['timestepping-timesteps.tex']);

h(end+1) = figure; hold on;
% Previously incorrectly used cumsum. This is already cumulative timesteps
% for each collocation point!!!!!
% n_steps_approx = cumsum(sum(n_steps_in_super.*I_in_super,2,'omitnan'));
n_steps_approx = (sum(n_steps_in_super.*I_in_super,2,'omitnan'));
% n_steps_estimation = cumsum(sum(n_steps_in_super.*I_star_in_super,2,'omitnan') - sum(n_steps_in_super.*I_in_super,2,'omitnan'));
n_steps_estimation = (sum(n_steps_in_super.*I_star_in_super,2,'omitnan') - sum(n_steps_in_super.*I_in_super,2,'omitnan'));
plot(data_table{:,'t'}, n_steps_approx);
plot(data_table{:,'t'}, n_steps_estimation);
%matlab2tikz('filename',['cumulative-steps.tex']);

h(end+1) = figure; hold on;
plot(data_table{inds_valid,'t'},MI_max(inds_valid,:));
set(gca,'XScale','log');
legend('show');
%matlab2tikz('filename',['max-mi.tex']);

% disp(MI{end});

savefig(h,'outputfigures.fig');


%%
pi_ts_simple = mean(ge_est_in_super(:,:).*I_in_super(:,:),2,'omitnan');
pi_ts_simple2 = ge_zero.';
pi_simple = data_table{:,'pi_I'} + pi_ts_simple;
pi_simple2 = data_table{:,'pi_I'} + pi_ts_simple2;

estsimple_union = interp1(data_table{inds_valid,'t'}, pi_simple(inds_valid), t_union,'linear');
estsimple2_union = interp1(data_table{inds_valid,'t'}, pi_simple2(inds_valid), t_union,'linear');
efficiency_full = est_union./error_union;
efficiency_simple = estsimple_union./error_union;
efficiency_simple2 = estsimple2_union./error_union;


%% Put blank lines in datatable
jj=1;
for ii = 1:size(data_table,1);
    data_table_blanks(jj,:) = data_table(ii,:);
    jj =jj +1;
    if inds_valid(ii) == 0
%         data_table_blank{jj,1} = ';
        jj=jj+1;
    end
end

%% Refinement times
data_t_shifted = [0 ;data_table{1:end,'t'}];
refinement_times = data_t_shifted(~inds_valid);

%% Refinement indices
refinement_table = data_table(~inds_valid,{'J','RM','pi_I_alpha'});
refinement_inds=[];
for ii = 1:size(refinement_table,1)
    refinement_inds_ii = refinement_table{ii,'J'}{1};
    refinement_inds = [refinement_inds; repmat(refinement_times(ii),[size(refinement_inds_ii,1),1]) refinement_inds_ii];
end

%% Refinement error estimators
refinement_estimators = [];
top10 = [];
for ii = 1:size(refinement_table,1)
    RM_ii = refinement_table{ii,'RM'}{1};
    pi_I_alpha_ii = refinement_table{ii,'pi_I_alpha'}{1};
    [~, sorted_inds] = sort(pi_I_alpha_ii,'descend');
    if length(sorted_inds) > 9
        top10_inds = sorted_inds(1:10);
    else
        top10_inds= sorted_inds;
    end
    pi_I_alpha_ii_sorted = pi_I_alpha_ii(sorted_inds);
    pi_I_alpha_10_ii = pi_I_alpha_ii(top10_inds);
    top10 = [top10; repmat(refinement_times(ii),[length(top10_inds),1])  RM_ii(top10_inds,:), pi_I_alpha_10_ii(:)];
    refinement_estimators = [refinement_estimators; repmat(refinement_times(ii),[size(RM_ii,1),1]) RM_ii(sorted_inds,:) pi_I_alpha_ii_sorted(:)];
end

%% nColloc
n_colloc_approx = sum(I_in_super,2,'omitnan');
n_colloc_estimation = sum(I_star_in_super,2,'omitnan') - n_colloc_approx;

%% Write matrices
writematrix([data_table{inds_valid,'t'}, mean(ge_est_in_super(inds_valid,:).*I_in_super(inds_valid,:),2,'omitnan'), min(ge_est_in_super(inds_valid,:).*I_in_super(inds_valid,:),[],2,'omitnan'), max(ge_est_in_super(inds_valid,:).*I_in_super(inds_valid,:),[],2,'omitnan'), ge_zero(inds_valid).'],'ge_est.dat','Delimiter','space');
writematrix([t_union, mean(ge_ref_in_super,2,'omitnan'), min(ge_ref_in_super,[],2,'omitnan'), max(ge_ref_in_super,[],2,'omitnan')],'ge_ref.dat','Delimiter','space');
writematrix([t_union, mean(efficiency_ge,2,'omitnan'), min(efficiency_ge,[],2,'omitnan'), max(efficiency_ge,[],2,'omitnan')],'ge_eff.dat','Delimiter','space');
writematrix([data_table{:,{'t','pi','pi_I','pi_I_delta','pi_delta','E'}}, pi_simple, pi_ts_simple, pi_simple2, pi_ts_simple2],'components-errorest.dat','Delimiter','space');
writematrix([data_table_blanks{:,{'t','pi','pi_I','pi_I_delta','pi_delta','E'}}],'components-errorest-blanks.dat','Delimiter','space'); 
writematrix([reference.data_table{:,{'t','error'}}, est_union, estsimple_union, estsimple2_union, efficiency_full, efficiency_simple, efficiency_simple2, error_union],'error.dat','Delimiter','space');
writematrix([data_table{:,{'t'}}, dt_mean.',dt_std.', dt_min.', dt_max.', dt_std.'./dt_mean.'], 'timesteps.dat','Delimiter','space');
writematrix([data_table{inds_valid,'t'},MI_max(inds_valid,:)], 'max-mi.dat','Delimiter','space');
writematrix([data_table{:,'t'},data_table{:,'delta_t'}],'alg-timesteps.dat','Delimiter','space') 
writematrix([data_table{:,'t'}, n_steps_approx, n_steps_estimation, n_steps_approx + n_steps_estimation, n_steps_zero(:)],'nsteps.dat','Delimiter','space');
writematrix([data_table{:,'t'}, n_colloc_approx, n_colloc_estimation, n_colloc_approx+n_colloc_estimation],'ncolloc.dat','Delimiter','space');
writematrix([refinement_times(:)],'refinement_times.dat','Delimiter','space');
writematrix([refinement_inds],'refinement_inds.dat','Delimiter','space');
writematrix([refinement_estimators],'refinement_estimators.dat','Delimiter','space');
writematrix([top10],'top10.dat','Delimiter','space');
end

function [e_mean,e_std, e_min, e_max] = plot_error_stats(axHandle, time, errors,colour)
for ii = 1:length(errors);
    [e_mean(ii), e_std(ii), e_min(ii), e_max(ii)] = extract_summary_stats(errors{ii});
end
plot(axHandle,time,e_mean,'-','Color',colour);
plot(axHandle,time,e_mean + e_std,'--','Color',colour);
plot(axHandle,time,e_mean - e_std,'--','Color',colour);
plot(axHandle,time,e_min,':','Color',colour);
plot(axHandle,time,e_max,':','Color',colour);
end