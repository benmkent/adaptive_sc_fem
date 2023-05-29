function [R] = filter_rows_columns(M, filter_rows,filter_columns)
R = M(filter_rows,filter_columns);
end