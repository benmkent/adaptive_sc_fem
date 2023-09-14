addpath ../../adaptive_sc_fem/;
problem = define_problem('doubleglazing-timebc');

%% Set up reference
params = define_params('smolyak-ref');
params.letol = 1e-5;
adaptive_sc_fem;
save(['reference.mat'],'reference','data_table','fem','problem','params', '-v7.3')

%% Run experiments
% names = {'smolyak-1','smolyak-2','smolyak-3','smolyak-4','smolyak-5'};
% names = {'smolyak-3','smolyak-4','smolyak-5'};

% problem = define_problem('doubleglazing-206');

reference = define_reference('test-folder');
params = define_params('adaptive-base');
params.letol = 1e-4;

%     params.plot = 1;

adaptive_sc_fem;
save(['adaptive.mat'],'reference','data_table','fem','problem','params','-v7.3')

%% Post process
% figOrig = figure(1); hold on;
% axOrig = gca(figOrig);
% set(axOrig,'XScale','log','YScale','log')
% figRes = figure(2); hold on;
% axRes = gca(figRes);
% set(axRes,'XScale','log','YScale','log')
% figErrorSum = figure(3); hold on;
% axErrorSum = gca(figErrorSum);
% set(axErrorSum,'XScale','log','YScale','log')
%
% c_array = lines(length(names));

% ii_c = 1;
% for name = names
%     c = c_array(ii_c,:);
%     load(name{1})
%     BK original estimation
%         plot(axOrig,reference.data_table{:,'t'},reference.data_table{:,'error'},'-x','Color',c);
%     plot(axOrig,data_table{:,'t'},data_table{:,'pi_delta'},'-o','Color',c);
%
%     Residual estimation
%         error_L2T_H=  sqrt(cumtrapz([0; reference.data_table{:,'t'}],[0; fem.lambda * reference.data_table{:,'error_H'}].^2));
%         error_X= sqrt(cummax([0;reference.data_table{:,'error'}]).^2 + cumtrapz([0; reference.data_table{:,'t'}],[0; fem.lambda * reference.data_table{:,'error_H'}].^2));
%
%         plot(axRes,[0;reference.data_table{:,'t'}], [error_X] ,'-x','Color',c);
%         plot(axRes,data_table{:,'t'},data_table{:,'pi_t'},'-o');
%         plot(axRes,[0;data_table{:,'t'}],[0;data_table{:,'pi_t'}],'-o','Color',c);
%
%         plot(axErrorSum, [0;reference.data_table{:,'t'}], [0;(reference.data_table{:,'error'})],'--o','Color',c);
%         plot(axErrorSum, [0;reference.data_table{:,'t'}], [error_L2T_H], '--x','Color',c);
%         plot(axErrorSum, [0;reference.data_table{:,'t'}], [(error_X)],'Color',c);
%
%         ii_c = ii_c + 1;
%
%     figure;
%     I_final = data_table{end,'I'};
%     idx_final = reshape([I_final{1}.idx],[],20);
%     [~,full_idx] = check_set_admissibility(idx_final);
%     Set up four large recirculating, then 16 small
%     sqrt_n = 2; %round(sqrt(n_dim));
%     dist_centres = 2/sqrt_n;
%     centres = -1 + dist_centres*0.5 + dist_centres * (0:(sqrt_n-1));
%     figure;
%     heatmap(centres,centres,reshape(max(full_idx(:,1:4)),[2,2]));
%     Set up four large recirculating, then 16 small
%     sqrt_n = 4; %round(sqrt(n_dim));
%     dist_centres = 2/sqrt_n;
%     centres = -1 + dist_centres*0.5 + dist_centres * (0:(sqrt_n-1));
%     figure;
%     heatmap(centres,centres,reshape(max(full_idx(:,5:20)),[4,4]));
%     set(gca,'xdir','normal');
%     set(gca,'ydir','normal');
% end