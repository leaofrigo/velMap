xxCentroLin = linspace(min(xxCentro),max(xxCentro),250);
cellDepthsLin = linspace(min(cellDepths),max(cellDepths),50);
[xxCentroGrid,cellDepthsGrid] = meshgrid(xxCentroLin,cellDepthsLin);

AAgridVelU = griddata(xxCentro,cellDepthsCentro,VelUCurl,xxCentroGrid,cellDepthsGrid);
figure
hp=pcolor(xxCentroGrid,cellDepthsGrid,AAgridVelU);
set(hp,'EdgeColor','none')
colorbar
hold on
plot(tracknew1, depthnew1)
xlim([tracknew1(1),tracknew1(end)])
ylim([min(depth)+.1*diff([0,min(depth)]),0]);


% imagesc(AAgrid)