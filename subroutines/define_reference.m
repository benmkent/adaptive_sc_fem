function reference = define_reference(varargin)
reference.times = logspace(-2,2,20);
reference.data_table = table([],[],[],[],[],[],[],[],'VariableNames',{'t','u','I','error','error_H','ts_error','exp','std'});
reference.save = 0;

if nargin == 0
    load_file = default('Load file reference.mat',1);
    if load_file == 1
        try
            reference_approximation = load('reference.mat');
        catch e
            error(e.message);
        end
        reference.reference_table = reference_approximation.reference.data_table;
        reference.times = reference_approximation.reference.times;
        reference.data_table = [];
        reference.fem = reference_approximation.fem;
    else
        reference.reference_table = [];
    end
else
    switch varargin{1}
        case 'test-folder'
            try
            reference_approximation = load('reference.mat');
        catch e
            error(e.message);
        end
        reference.reference_table = reference_approximation.reference.data_table;
        reference.times = reference_approximation.reference.times;
        reference.fem = reference_approximation.fem;
    end
end