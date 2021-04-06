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
    AAgridVelU(TFTemp2) = nan;    
end

% [PerimBotX,PerimBotY] = meshgrid(PerimBotX,PerimBotY);
% [PerimSupX,PerimSupY] = meshgrid(PerimSupX,PerimSupY);


