%figure; 
plot(reference.data_table{:,'t'},reference.data_table{:,'error'},'DisplayName','ptwise')
hold on; plot([0;reference.data_table{:,'t'}],sqrt(cumtrapz([0;reference.data_table{:,'t'}],[0;reference.data_table{:,'error'}].^2)),'DisplayName','L^2')
plot(data_table{:,'t'}, data_table{:,'pi_x'},'DisplayName','x','Marker','+')
plot(data_table{:,'t'}, data_table{:,'pi_t'},'DisplayName','t','Marker','o')
plot(data_table{:,'t'}, data_table{:,'pi_y'},'DisplayName','y','Marker','x')
set(gca(),'YScale','log')
set(gca(),'XScale','log')
legend('show')