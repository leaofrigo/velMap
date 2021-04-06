Space = 100;
xmin = min(Quiver3VectorXTot(:)); 
ymin = min(Quiver3VectorYTot(:)); 
zmin = min(Quiver3VectorZTot(:));

xmax = max(Quiver3VectorXTot(:)); 
ymax = max(Quiver3VectorYTot(:)); 
zmax = max(Quiver3VectorZTot(:));

x1 = linspace(xmin,xmax,Space);
y1 = linspace(ymin,ymax,Space);
z1 = linspace(zmin,zmax,Space);
[X1,Y1,Z1] = meshgrid(x1,y1,z1);
Quiver3VectorMagTot = sqrt(Quiver3VectorUTot.^2 +Quiver3VectorVTot.^2+Quiver3VectorWTot.^2 );
VEL3D = griddata(Quiver3VectorXTot,Quiver3VectorYTot,Quiver3VectorZTot,Quiver3VectorMagTot,X1,Y1,Z1);
% YrRR = YrR;
% YrRR(isnan(FrR)) = nan;
% nanmean(YrRR),XrR(1,:);
h1 = figure;
Gr = surf(YrR,XrR,FrR-A.GPS.Altitude(1));
set(Gr,'EdgeColor', 'none');
Xlim = get(gca,'Xlim'); Ylim =get(gca,'Ylim');Zlim = get(gca,'Zlim');
i=1;
j=1;
[LowerLim, UpperLim,s] = limites(VEL3D);
ylabel('Latitude (deg)')
xlabel('Longitude (deg)')
zlabel('Profundidade (m)')
% h2 = figure;
% Grad = gradient(nanmean(YrRR),XrR(1,:));
% xd = zeros(Space,Space,Space);
% yd=xd;zd=xd;
% Zmesh = linspace(min(FrR(:))-A.GPS.Altitude(1),max(FrR(:))-A.GPS.Altitude(1),Space);
% Xmesh = linspace(min(YrR(:)),max(YrR(:)),Space);
% meshX = meshgrid(Xmesh,Xmesh);
% meshZ = meshgrid(Zmesh,Zmesh);
% for i = 1:Space
%     clf
%     hslice = surf(meshX,ones(100)*XrR(1,i),meshZ');
%     rotate(hslice,[0,0,1],atan(Grad(i)))
%     xd(:,:,i) = get(hslice,'XData');
%     yd(:,:,i) = get(hslice,'YData');
%     zd(:,:,i) = get(hslice,'ZData');
% 
%     
% end
% close(figure(h2))
% i=1;
% figure(h1)
while 1
    
    hold on
%     fe = slice(X1,Y1,Z1,VEL3D,xd(:,:,i),yd(:,:,i),zd(:,:,i));
    fe = slice(X1,Y1,Z1,VEL3D,0,XrR(1,i),0);
    if j == 1
        colorbar('eastoutside')
        j = 0;
        h_bar = findobj(gcf,'Tag','Colorbar');
        set(get(h_bar,'xlabel'),'String', 'Velocidade(m/s)');
    end
    caxis([LowerLim, UpperLim]);
    set(fe,'EdgeColor', 'none');
    axis([Xlim Ylim Zlim]);
    i=i+1;
    if i>Space
        i=1;
    end    
    pause(0.1);
    delete(fe)
    
    
end

