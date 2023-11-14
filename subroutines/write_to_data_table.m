function data_table = write_to_data_table(params,varargin)
%WRITE_TO_DATA_TABLE Writes data structures to data table
% Inputs
%    params     initialises data_table if only argument
%    varargin   datatlable and the data to be written to the table.
% output
%    data_table initialised or updated data_table

if nargin == 1
    switch params.adapt_type
        case 'hierarchical'
            data_table = table([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],...
                'VariableNames',{...
                't', ...
                'delta_t', ...
                'E', ...
                'pi', ...
                'pi_I', ...
                'pi_I_delta', ...
                'pi_delta', ...
                'pi_I_alpha', ...
                'u_z',...
                'RM', ...
                'J', ...
                'ge_estimate', ...
                'dt_z', ...
                'n_steps', ...
                'n_steps_lofi', ...
                'I_star', ...
                'I'});
        case 'residual'
            data_table = table([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],...
                'VariableNames',{...
                't', ...
                'delta_t', ...
                'E', ...
                'pi_x',...
                'pi_t',...
                'pi_y',...
                'pi_x_r',...
                'pi_x_z_T',...
                'u_z',...
                'pi_t_r',...
                'eta_t_z_du',...
                'eta_t_z_u_w',...
                'eta_t_z_eta_w',...
                'pi_y_r',...
                'pi_xty',...
                'pi_y_mi_r',...
                'y_mi_r',...
                'RM', ...
                'J', ...
                'dt_z', ...
                'n_steps', ...
                'I_star', ...
                'I'});
    end
else
    data_table = varargin{1};
    t = varargin{2};
    delta_t= varargin{3};
    E= varargin{4};
    error= varargin{5};
    RM= varargin{6};
    J= varargin{7};
    ge_estimates= varargin{8};
    dt_z_tplusdt= varargin{9};
    n_steps= varargin{10};
    n_steps_lofi= varargin{11};
    I_star = varargin{12};
    I = varargin{13};
    switch params.adapt_type
        case 'hierarchical'

            data_table = [data_table;
                table([t], ...
                [delta_t], ...
                [E], ...
                [error.pi], ...
                [error.pi_I], ...
                [error.pi_I_delta], ...
                [error.pi_delta], ...
                {error.pi_I_alpha}, ...
                {RM}, ...
                {J}, ...
                {ge_estimates}, ...
                {dt_z_tplusdt}, ...
                {n_steps}, ...
                {n_steps_lofi}, ...
                {I_star}, ...
                {I},...
                'VariableNames',data_table.Properties.VariableNames)...
                ];
        case 'residual'

            data_table = [data_table;
                table([t], ...
                [delta_t], ...
                [E], ...
                error.pi_x,...
                error.pi_t,...
                error.pi_y,...
                error.pi_x_r,...
                error.eta_x_z_T,...
                error.eta_u_z,...
                error.pi_t_r,...
                error.eta_t_z_du,...
                error.eta_t_z_u_w,...
                error.eta_t_z_eta_w,...
                error.pi_y_r,...
                error.pi_xty,...
                {error.pi_y_mi_r},...
                {error.y_mi_r},...
                {RM}, ...
                {J}, ...
                {dt_z_tplusdt}, ...
                {n_steps}, ...
                {I_star}, ...
                {I},...
                'VariableNames',data_table.Properties.VariableNames)...
                ];

    end
end
