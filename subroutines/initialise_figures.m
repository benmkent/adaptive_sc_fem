function axArray = initialise_figures()
figExp = figure(1); axArray.axExp = gca(figExp);
figStd = figure(2); axArray.axStd = gca(figStd);
figError = figure(3); axArray.axError = gca(figError);
figSnap = figure(4); axArray.axSnap = gca(figSnap);
axArray.figSnap = figSnap;
axArray.axSnap1 = subplot(2,2,1);
axArray.axSnap2 = subplot(2,2,2);
axArray.axSnap3 = subplot(2,2,3);
axArray.axSnap4 = subplot(2,2,4);
figVid = figure(5); axArray.axVid = gca(figVid);
figEtaT = figure(6); axArray.axEtaT = gca(figEtaT);
figTimesteps = figure(7); axArray.axTimesteps = gca(figTimesteps);
figMaxLevel = figure(8); axArray.axMaxLevel = gca(figMaxLevel);
end