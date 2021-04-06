numofvalues = 100;
x = linspace(min(GridVar.Longitudeidx),max(GridVar.Longitudeidx),numofvalues);
y = linspace(min(GridVar.Latitudeidx),max(GridVar.Latitudeidx),numofvalues);
[X,Y] = meshgrid(x,y);
Ve = griddata(GridVar.Longitudeidx,GridVar.Latitudeidx,GridVar.VelocityEidx,X,Y);
Vn = griddata(GridVar.Longitudeidx,GridVar.Latitudeidx,GridVar.VelocityNidx,X,Y);
quiver(X,Y,Ve,Vn)




h = streamline(X,Y,Ve,Vn,-49.3355,-4.9764);