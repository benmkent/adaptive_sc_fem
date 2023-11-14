function  output_table = post_process(data_table, t, delta_t, Z_I_star, Z_I_star_lofi,...
    I_star, I, E, error, work, RM, J, problem, params, fem, precompute, ...
    axArray)
%POST_PROCESS   Post process approximation errors to data table and plot
%
%   Inputs      data_table     input data_table
%               t              t at start of interval
%               delta_t        algorithm timestep
%               Z_I_star       collocation point data structure
%               Z_I_star_lofi  lofi collocation point data structure
%               I_star          Enhanced sparse grid
%               I               Sparse grid
%               E               Error tolerance
%               error           error indicator structure
%               work            Computed work for each MI
%               RM              Reduced margin (or margin in residual case)
%               J               Marked miset
%               problem         problem structure
%               params          parameter structure
%               fem             Finite element structure
%               precompute      Precomputed data
%               axArray         Figure axes handles

%% Set up
I_star_r = reduce_sparse_grid(I_star);
I_r = reduce_sparse_grid(I);
if I_star_r.size == I_r.size
    both_I_star = 1:I_r.size;
else
    [I_star_only, both_I_star, both_I, I_only] = compare_sparse_grids(I_star, I_star_r, I, I_r);
end
%% Extract data
ge_estimates = zeros(length(Z_I_star),1);
dt_z_tplusdt = zeros(length(Z_I_star),1);
for ii = 1:length(Z_I_star)
    ge_estimates(ii) = Z_I_star{(ii)}.ge_estimate;
    dt_z_tplusdt(ii) = Z_I_star{(ii)}.dt_z_tplusdt;
    n_steps(ii) = Z_I_star{(ii)}.n_steps;
    switch params.adapt_type
        case 'hierarchical'
            n_steps_lofi(ii) = Z_I_star_lofi{(ii)}.n_steps;
        case 'residual'
            n_steps_lofi(ii) = 0;
    end
end

%% Write to data table
data_table = write_to_data_table(params,data_table, t, delta_t, E, error, RM, J, ge_estimates, dt_z_tplusdt, n_steps, n_steps_lofi, I_star, I);

%% Plot data
if params.plot == 1
    axExp = axArray.axExp;
    axStd = axArray.axStd;
    axError = axArray.axError;
    figSnap = axArray.figSnap;
    axSnap1 = axArray.axSnap1;
    axSnap2 = axArray.axSnap2;
    axSnap3 = axArray.axSnap3;
    axSnap4 = axArray.axSnap4;
    axEtaT = axArray.axEtaT;
    axVid = axArray.axVid;
    axTimesteps = axArray.axTimesteps;
    axMaxLevel = axArray.axMaxLevel;

    u_z_tplusdt= zeros(size(Z_I_star{1}.u_z_tplusdt,1),length(both_I_star));
    for ii = 1:length(both_I_star)
        u_z_tplusdt(:,ii) = abs(Z_I_star{both_I_star(ii)}.u_z_tplusdt);
        z_tplusdt(:,ii) = Z_I_star{both_I_star(ii)}.z;
    end
    vec_exp_u = sum([u_z_tplusdt] .* I_r.weights,2);
    plot_surf(axExp,vec_exp_u, fem, params);
    caxis(axExp,[0,1]);
    title(axExp,['Exp t=' num2str(t)]);

    u_z_interpolated = interpolate_on_sparse_grid(I,I_r,u_z_tplusdt, I_star_r.knots);
    vec_exp_u_sq = sum(abs(u_z_interpolated).^2 .* I_star_r.weights,2);
    vec_var = vec_exp_u_sq - vec_exp_u.^2;
    vec_var(vec_var < 0) = 0;
    plot_surf(axStd,real(sqrt(vec_var)),fem,params);
    caxis(axStd,[0,0.1]);
    %     caxis(axStd,[-9,0]);
    title(axStd,['Std t=' num2str(t)]);

    set(axTimesteps,'NextPlot','add');
    errorbar(axTimesteps,t, mean(dt_z_tplusdt),std(dt_z_tplusdt));
    scatter(axTimesteps,t, length(Z_I_star));
    set(axTimesteps,'XScale','log','YScale','log');
    xlabel(axTimesteps,'t')
    legend(axTimesteps,'dt','nColloc')


    ind = round(linspace(1,size(u_z_tplusdt,2),4));
    mi = get_mi_set(I);
    max_inds = max(mi,[],1);
    %
    %     [~,ind_max_level] = max(max_inds);
    %
    %     [row,col] = find(mi == max_inds);
%     if isempty(J)
        times_to_plot = params.plot_times(params.plot_times >= t-delta_t & params.plot_times < t);
%         data_1 = interp1(Z_I_star{ind(1)}.t_z,Z_I_star{ind(1)}.u_z.',times_to_plot(:));
%         data_2 = interp1(Z_I_star{ind(2)}.t_z,Z_I_star{ind(2)}.u_z.',times_to_plot(:));
%         data_3 = interp1(Z_I_star{ind(3)}.t_z,Z_I_star{ind(3)}.u_z.',times_to_plot(:));
%         data_4 = interp1(Z_I_star{ind(4)}.t_z,Z_I_star{ind(4)}.u_z.',times_to_plot(:));
%         writematrix([times_to_plot(:),data_1],'data1.dat','WriteMode','append');
%         writematrix([times_to_plot(:),data_2],'data2.dat','WriteMode','append');
%         writematrix([times_to_plot(:),data_3],'data3.dat','WriteMode','append');
%         writematrix([times_to_plot(:),data_4],'data4.dat','WriteMode','append');
%     end
%     for ii_t = length(times_to_plot)
        plot_surf(axSnap1,u_z_tplusdt(:,ind(1)),fem,params)
        caxis(axSnap1,[0,1]);
        title(axSnap1,['y=' num2str(z_tplusdt(1,ind(1))) ', t=' num2str(t)]);

        plot_surf(axSnap2,u_z_tplusdt(:,ind(2)),fem,params)
        caxis(axSnap2,[0,1]);
        title(axSnap2,['y=' num2str(z_tplusdt(1,ind(2))) ', t=' num2str(t)]);

        plot_surf(axSnap3,u_z_tplusdt(:,ind(3)),fem,params)
        caxis(axSnap3,[0,1]);
        title(axSnap3,['y=' num2str(z_tplusdt(1,ind(3))) ', t=' num2str(t)]);

        plot_surf(axSnap4,u_z_tplusdt(:,ind(4)),fem,params)
        caxis(axSnap4,[0,1]);
        title(axSnap4,['y=' num2str(z_tplusdt(1,ind(4))) ', t=' num2str(t)]);
%     end

    switch params.adapt_type
        case 'hierarchical'
            set(axError,'NextPlot','replaceall');
            plot(axError,data_table,'t','E')
            set(axError,'NextPlot','add');
            plot(axError,data_table,'t','pi')
            plot(axError,data_table,'t','pi_I')
            plot(axError,data_table,'t','pi_I_delta')
            plot(axError,data_table,'t','pi_delta')
            set(axError,'XScale','log','YScale','log');
            legend(axError,'show','Location','southeast');
            hold off;

            pi_alpha = error.pi_I_alpha;
            if ~isempty(pi_alpha)
                cla(axMaxLevel);
                bar(axMaxLevel,max_inds);
                plot(axMaxLevel,pi_alpha);
                drawnow();
            end
        case 'residual'
            set(axEtaT,'NextPlot','replaceall');
            plot(axEtaT,data_table,'t','E')
            set(axEtaT,'NextPlot','add');
            plot(axEtaT,data_table,'t','pi_x_r')
            plot(axEtaT,data_table,'t','pi_t_r')
            plot(axEtaT,data_table,'t','pi_y_r')
            set(axEtaT,'XScale','log','YScale','log');
            legend(axEtaT,'show','Location','southeast');
            hold off;
            set(axError,'NextPlot','replaceall');
            set(axError,'NextPlot','add');
            plot(axError,data_table,'t','pi_x')
            plot(axError,data_table,'t','pi_t')
            plot(axError,data_table,'t','pi_y')
            plot(axError,data_table,'t','pi_xty')
            set(axError,'XScale','log','YScale','log');
            legend(axError,'show','Location','southeast');
            hold off;
    end

    drawnow();

end

%% Return data_table
output_table = data_table;