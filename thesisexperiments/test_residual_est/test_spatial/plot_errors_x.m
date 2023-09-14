figure; 
amin=0.1;

load test-spatial-l2.mat
plot(reference.data_table{:,'t'},reference.data_table{:,'error'},'DisplayName','ptwise')
hold on; plot([0;reference.data_table{:,'t'}],sqrt(amin*cumtrapz([0;reference.data_table{:,'t'}],[0;reference.data_table{:,'error_H'}].^2)),'DisplayName','L^2')
plot(data_table{:,'t'}, data_table{:,'pi_x'},'-x','DisplayName','x')
plot(data_table{:,'t'}, data_table{:,'pi_t'},'-o','DisplayName','t')
plot(data_table{:,'t'}, data_table{:,'pi_y'},'-+','DisplayName','y')
set(gca(),'YScale','log')
set(gca(),'XScale','log')
legend('show')
load test-spatial-l3.mat
plot(reference.data_table{:,'t'},reference.data_table{:,'error'},'DisplayName','ptwise')
hold on; plot([0;reference.data_table{:,'t'}],sqrt(amin*cumtrapz([0;reference.data_table{:,'t'}],[0;reference.data_table{:,'error_H'}].^2)),'DisplayName','L^2')
plot(data_table{:,'t'}, data_table{:,'pi_x'},'-x','DisplayName','x')
plot(data_table{:,'t'}, data_table{:,'pi_t'},'-o','DisplayName','t')
plot(data_table{:,'t'}, data_table{:,'pi_y'},'-+','DisplayName','y')
set(gca(),'YScale','log')
set(gca(),'XScale','log')
legend('show')
load test-spatial-l4.mat
plot(reference.data_table{:,'t'},reference.data_table{:,'error'},'DisplayName','ptwise')
hold on; plot([0;reference.data_table{:,'t'}],sqrt(amin*cumtrapz([0;reference.data_table{:,'t'}],[0;reference.data_table{:,'error_H'}].^2)),'DisplayName','L^2')
plot(data_table{:,'t'}, data_table{:,'pi_x'},'-x','DisplayName','x')
plot(data_table{:,'t'}, data_table{:,'pi_t'},'-o','DisplayName','t')
plot(data_table{:,'t'}, data_table{:,'pi_y'},'-+','DisplayName','y')
set(gca(),'YScale','log')
set(gca(),'XScale','log')
legend('show')
% 
% figure;
% plot_data_tifiss(1,1,data_table{20,'u_z'}{1},data_table{20,'pi_x_z_T'}{1},fem.ev,fem.xy)