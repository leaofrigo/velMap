function p = SmoothFull(xxCentro,NumOfCellsX,cellDepths,cellDepthsCentro,NumOfCellsY,VelVector,...
    Sumnum2,xx,Direcao,C,StartEdge,tracknew1,depthnew1,OpenFiles,XXQuiver,YYQuiver,...
    Quiver3VectorU,Quiver3VectorV,Quiver3VectorW,Quiver3VectorMag,Quiver3VectorErr,VelRL,Directory,...
    quiv,handles)

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
clf
p=pcolor(xxCentroGrid,cellDepthsGrid,AAgridVel);
set(p,'EdgeColor','none')
colorbar('south')
h_bar = findobj(gcf,'Tag','Colorbar');
set(get(h_bar,'xlabel'),'String', 'Velocity(m/s)');
initpos = get(h_bar,'Position');
initfontsize = get(h_bar,'FontSize');
set(h_bar,'Position',[initpos(1)+initpos(3)/2 initpos(2) initpos(3)/2 initpos(4)/2],...
    'FontSize',initfontsize*.75)
hold on

if handles.run3d == 1
    Pos = get(handles.Options3D.plotpostion3d,'UserData');
else
    Pos = get(handles.Options2D.plotpostion,'UserData');
end
set(gcf, 'Position', Pos);
plot(tracknew1, depthnew1)
colormap(flipud(C([1:4:16,14:50,50:2:end-125,...
    end-125:end-50,end-50:4:end-30],:))) % make the color bar look like the River Servoyor Live
xlim([tracknew1(1),tracknew1(end)])
ylim([min(depthnew1)+.1*diff([0,min(depthnew1)]),0]);
xlabel('Length (m)')
ylabel('Depth (m)')
if StartEdge == 1
    set(gca,'XDir','reverse')
end
% if n > 1
%     OpenFiles = OpenFiles{n};
% end
if quiv{1} ==1
    hold on
   switch quiv{2}
        case 'North'
            Vel1 = Quiver3VectorV;      
        case 'East'
            Vel1 = Quiver3VectorU;        
        case 'Upstream'
            Vel1 = Quiver3VectorW;        
        case 'Magnitude'
            Vel1 = Quiver3VectorMag;
        case 'Error'
            Vel1 = Quiver3VectorErr; 
        case 'Transverse'
            Vel1 = VelRL(:,1);
        case 'Longitudinal'     
            Vel1 = VelRL(:,2);
    end
    switch quiv{3}
        case 'North'
            Vel2 = Quiver3VectorV;      
        case 'East'
            Vel2 = Quiver3VectorU;        
        case 'Upstream'
            Vel2 = Quiver3VectorW;        
        case 'Magnitude'
            Vel2 = Quiver3VectorMag;
        case 'Error'
            Vel2 = Quiver3VectorErr; 
        case 'Transverse'
            Vel2 = VelRL(:,1);
        case 'Longitudinal'     
            Vel2 = VelRL(:,2);
    end
    hold on
    if handles.run3d == 1;
        quivScale = get(handles.Options3D.quiv3dOption,'UserData');
    else
        quivScale = get(handles.Options2D.quiv2dOption,'UserData');
    end
    if StartEdge ==1
        quiver(XXQuiver,YYQuiver,-Vel1,Vel2,quivScale)
    else
        quiver(XXQuiver,YYQuiver,Vel1,Vel2,quivScale)
    end
end
    
switch Direcao
    case 'North'
        title('Velocity North');
        [LowerVN,UpperVN] = limites(VelVector); 
        caxis([LowerVN UpperVN]);
        savefig(gcf,[Directory 'Section ' OpenFiles '\Velocity_North_Smooth.fig']);
    case 'East'
        title('Velocity East');
        [LowerVE,UpperVE] = limites(VelVector); 
        caxis([LowerVE UpperVE]);
        savefig(gcf,[Directory 'Section ' OpenFiles '/Velocity_East_Smooth.fig']);
    case 'Upstream'
        title('Velocity Upstream');
        [LowerVU,UpperVU] = limites(VelVector); 
        caxis([LowerVU UpperVU]);
        savefig(gcf,[Directory 'Section ' OpenFiles '/Velocity_Upstream_Smooth.fig']);
    case 'Error'
        title('Error Velocity');
        [LowerVD,UpperVD] = limites(VelVector); 
        caxis([LowerVD UpperVD]);
        savefig(gcf,[Directory 'Section ' OpenFiles '/Velocity_Upstream_Smooth.fig']);
    case 'Magnitude'
        title('Velocity magnitude');
        [~,UpperVMag,StdMag] = limites(VelVector); 
        UpperVMag1 = ceil(max(max(VelVector))*10)/10;
        hold on
        savefig(gcf,[Directory 'Section ' OpenFiles '/Velocity_Magnitude_Smooth.fig']);
    case 'Transverse'
        [LowerVR,UpperVR] = limites(VelVector); 
        caxis([LowerVR UpperVR]);
        title('Velocity transverse')
        savefig(gcf,[Directory 'Section ' OpenFiles '/Velocity_Transversal_Smooth.fig']);
    case 'Longitudinal'
        title('Velocity Longitudinal')
        [LowerVL,UpperVL] = limites(VelVector); 
        caxis([LowerVL UpperVL]);
        savefig(gcf,[Directory 'Section ' OpenFiles '/Velocity_Longitudinal_Smooth.fig']);
end
