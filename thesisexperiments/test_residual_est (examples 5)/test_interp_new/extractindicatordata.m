amin=0.1;
plot(reference.data_table{:,'t'},reference.data_table{:,'error'},'DisplayName','ptwise')
hold on; plot([0;reference.data_table{:,'t'}],sqrt(amin*cumtrapz([0;reference.data_table{:,'t'}],[0;reference.data_table{:,'error_H'}].^2)),'DisplayName','L^2')
plot(data_table{:,'t'}, data_table{:,'pi_x'},'-x','DisplayName','x')
plot(data_table{:,'t'}, data_table{:,'pi_t'},'-o','DisplayName','t')
plot(data_table{:,'t'}, data_table{:,'pi_y'},'-+','DisplayName','y')

writecell([{'t','x','t','y','dx','dt','dy'}; num2cell([data_table{:,'t'},data_table{:,'pi_x'},data_table{:,'pi_t'},data_table{:,'pi_y'},...
    data_table{:,'pi_x'}./data_table{:,'delta_t'},data_table{:,'pi_t'}./data_table{:,'delta_t'},data_table{:,'pi_y'}./data_table{:,'delta_t'}])],...
    ['indicator-data-grid-' num2str(params.l_initial) '.dat'],'Delimiter','space');
writecell([{'t','L2','H1','E'};num2cell([[0;reference.data_table{:,'t'}],...
    [0;reference.data_table{:,'error'}],...
    sqrt(amin*cumtrapz([0;reference.data_table{:,'t'}],[0;reference.data_table{:,'error_H'}].^2)),...
    sqrt([0;reference.data_table{:,'error'}].^2 + amin*cumtrapz([0;reference.data_table{:,'t'}],[0;reference.data_table{:,'error_H'}].^2))])],...
    ['error-data-grid-' num2str(params.l_initial) '.dat'],'Delimiter','space');


nMargin = size(data_table{:,'y_mi_r'}{1},1);
nRows = size(data_table,1);
data_mi = reshape([data_table{:,'pi_y_mi_r'}{:}].',[nMargin,nRows]).';

writematrix([data_table{:,'t'},data_mi],['mi-indicators',num2str(params.l_initial),'.dat'],'Delimiter','space');
writematrix([data_table{:,'t'},data_mi./data_table{:,'delta_t'}],['mi-indicators-dt',num2str(params.l_initial),'.dat'],'Delimiter','space');
writematrix(data_table{:,'y_mi_r'}{1},['mi',num2str(params.l_initial),'.dat'],'Delimiter','space');