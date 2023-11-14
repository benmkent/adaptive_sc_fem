inputData = csvread('data1.dat');
vidObj = VideoWriter('data1.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 10;
open(vidObj);

figure;
ax = gca();
for ii = 1:size(inputData,1);
    plot_surf(ax,inputData(ii,2:end).',fem,params);
    zlim([-1e-1,1.5]);
    caxis([0,1.5]);
    title(num2str(inputData(ii,1)))
    frame = getframe(gcf);
    writeVideo(vidObj,frame)
end;