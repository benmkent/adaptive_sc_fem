figure(1)
% indices in 17 pt grid
i1_17=1;
i2_17=2;
i3_17=6;
i4_17=8;

I17 = reference.data_table{t1,'I'}{1};
I17_r = reduce_sparse_grid(reference.data_table{t1,'I'}{1});

t1=6;
i1 = i1_17;
i2=i2_17;
i3=i3_17;
i4=i4_17;

u01 = reference.data_table{t1,'u'}{1};
subplot(4,4,1);
contourf(reshape(u01(:,i1),[65,65]),40)
subplot(4,4,2);
contourf(reshape(u01(:,i2),[65,65]),40)
subplot(4,4,3);
contourf(reshape(u01(:,i3),[65,65]),40)
subplot(4,4,4);
contourf(reshape(u01(:,i4),[65,65]),40)

t2=11;
[~,old_inds_in_new] = compare_sparse_grids(reference.data_table{t2,'I'}{1},...
    reduce_sparse_grid(reference.data_table{t2,'I'}{1}),...
    I17,I17_r);
i1 = old_inds_in_new(i1_17);
i2 = old_inds_in_new(i2_17);
i3 = old_inds_in_new(i3_17);
i4 = old_inds_in_new(i4_17);
u1 = reference.data_table{t2,'u'}{1};
subplot(4,4,5);
contourf(reshape(u1(:,i1),[65,65]),40)
subplot(4,4,6);
contourf(reshape(u1(:,i2),[65,65]),40)
subplot(4,4,7);
contourf(reshape(u1(:,i3),[65,65]),40)
subplot(4,4,8);
contourf(reshape(u1(:,i4),[65,65]),40)

t3=14;
[~,old_inds_in_new] = compare_sparse_grids(reference.data_table{t3,'I'}{1},...
    reduce_sparse_grid(reference.data_table{t3,'I'}{1}),...
    I17,I17_r);
i1 = old_inds_in_new(i1_17);
i2 = old_inds_in_new(i2_17);
i3 = old_inds_in_new(i3_17);
i4 = old_inds_in_new(i4_17);
u10 = reference.data_table{t3,'u'}{1};
subplot(4,4,9);
contourf(reshape(u10(:,i1),[65,65]),40)
subplot(4,4,10);
contourf(reshape(u10(:,i2),[65,65]),40)
subplot(4,4,11);
contourf(reshape(u10(:,i3),[65,65]),40)
subplot(4,4,12);
contourf(reshape(u10(:,i4),[65,65]),40)

t4=20;
[~,old_inds_in_new] = compare_sparse_grids(reference.data_table{t4,'I'}{1},...
    reduce_sparse_grid(reference.data_table{t4,'I'}{1}),...
    I17,I17_r);
i1 = old_inds_in_new(i1_17);
i2 = old_inds_in_new(i2_17);
i3 = old_inds_in_new(i3_17);
i4 = old_inds_in_new(i4_17);
u100 = reference.data_table{t4,'u'}{1};
subplot(4,4,13);
contourf(reshape(u100(:,i1),[65,65]),40)
subplot(4,4,14);
contourf(reshape(u100(:,i2),[65,65]),40)
subplot(4,4,15);
contourf(reshape(u100(:,i3),[65,65]),40)
subplot(4,4,16);
contourf(reshape(u100(:,i4),[65,65]),40)
