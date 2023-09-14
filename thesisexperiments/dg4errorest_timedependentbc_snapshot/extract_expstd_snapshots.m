expmatrix = [reference.data_table.exp{:}];
stdmatrix = [reference.data_table.std{:}];
tmatrix = [reference.data_table.t];

writematrix([tmatrix(:), expmatrix.'],'expsnapshot.dat');
writematrix([tmatrix(:), stdmatrix.'],'stdsnapshot.dat');