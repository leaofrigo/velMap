function output = data2grid4multisess(xxCentro,NumOfCellsX,cellDepths,cellDepthsCentro,NumOfCellsY,VelVector,...
    Sumnum2,xx,tracknew1,depthnew1)

xxCentroLin = linspace(min(xxCentro(:)),max(xxCentro(:)),NumOfCellsX);
cellDepthsLin = linspace(min(cellDepths),max(cellDepths),NumOfCellsY);
[xxCentroGrid,cellDepthsGrid] = meshgrid(xxCentroLin,cellDepthsLin);

AAgridVel = griddata(xxCentro,cellDepthsCentro,VelVector,xxCentroGrid,cellDepthsGrid);

Yy = cellDepths(Sumnum2(1:length(Sumnum2)-1));
PerimSupX = zeros(1,2*(length(xx)-1));
PerimSupY = zeros(1,2*(length(xx)-1));
PerimSupY(1:2:end-1) = Yy;
PerimSupY(2:2:end) = Yy;
PerimSupX(1) = xx(1);
PerimSupX(end) = xx(end);
PerimSupX(2:2:end-2) = xx(2:end-1);
PerimSupX(3:2:end-1) = xx(2:end-1);

Yybottom = cellDepths(Sumnum2(2:end)-1);
PerimBotY = zeros(1,2*(length(xx)-1));
PerimBotY(1:2:end-1) = Yybottom;
PerimBotY(2:2:end) = Yybottom;
PerimBotX = PerimSupX;

for ij = 1:length(PerimBotX)-1    
    TFTemp = and(xxCentroGrid(1,:)>=PerimBotX(ij),xxCentroGrid(1,:)<=PerimBotX(ij+1));
    TFTemp1 = and(cellDepthsGrid(:,TFTemp)>=PerimBotY(ij),cellDepthsGrid(:,TFTemp)<=PerimSupY(ij));
    TFTemp2 = zeros(size(cellDepthsGrid));
    TFTemp2(:,TFTemp) = ~TFTemp1;
    TFTemp2 = logical(TFTemp2);
    cellDepthsGrid(TFTemp2) = nan; 
    xxCentroGrid(TFTemp2) = nan;
    AAgridVel(TFTemp2) = nan;    
end

output.xxCentro = xxCentroGrid;
output.cellDepthsGrid = cellDepthsGrid;
output.AAgridVel = AAgridVel;
output.track = tracknew1;
output.depth = depthnew1;
mindepth = min(depthnew1);
output.depthmin = mindepth;
