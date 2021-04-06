xGr = linspace(min(Quiver3VectorX),max(Quiver3VectorX),100);
yGr = linspace(min(Quiver3VectorY),max(Quiver3VectorY),100);
[XGR, YGR] = meshgrid(xGr,yGr);
GR = griddata(Quiver3VectorX,Quiver3VectorY,Quiver3VectorU,XGR, YGR);
pcolor(XGR,YGR,GR);
GRS = smooth2a(GR,50);
XGRS = smooth2a(XGR,50);
YGRS = smooth2a(YGR,50);
figure
pcolor(XGRS,YGRS,GRS);
% GRz = zeros(size(GRXY));
% slice(X1,Y1,Z1,VEL3D,[],y1(i),[]);