load('rf-construction.mat');
jj=3;
for ii=1:jj^2;
subplot(jj,jj,ii);
% contourf(xp,yp,-(1-xp.^2).*(1-yp.^2) + reshape(Vl(:,:)*2*(rand(64,1)-0.5),[129,129])),
contourf(xp,yp,-(1-xp.^2).*(1-yp.^2) + reshape(Vl(:,ii),[129,129])),

% caxis([-2,2]),
colorbar,
colormap('jet'),
axis square;
end
%%
figure
jj=8; inds=[1,2,3,4,40,50,60,64];
t = tiledlayout(2,4,'TileSpacing','Compact','Padding','Compact');
for ii=1:jj;
    nexttile
% ax = subplot(1,jj,ii);
% contourf(xp,yp,-(1-xp.^2).*(1-yp.^2) + reshape(Vl(:,:)*2*(rand(64,1)-0.5),[129,129])),
contourf(xp,yp,reshape(Vl(:,inds(ii)),[129,129]),5),
title(['i = ' num2str(inds(ii))]);
colorbar,
caxis([-0.5,0.5])
colormap('jet'),
axis square;

end