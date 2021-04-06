function p = CriarESalvarFigura(Face,VerticeXXYY,VelVector,tracknew1,...
    depthnew1,Direcao,Directory,OpenFiles,XXQuiver,YYQuiver,Quiver3VectorU,Quiver3VectorV,...
    Quiver3VectorW,Quiver3VectorMag,Quiver3VectorErr,VelRL,C,StartEdge,quiv,handles)
% if strcmp(Direcao,'Longitugional')==1
%     subplot(1,10,1:9)
% end
hold on
p = patch('Faces',Face,'Vertices',VerticeXXYY);
colorbar('south')
set(p,'FaceColor','flat',...
'FaceVertexCData',VelVector,...
'CDataMapping','scaled','EdgeColor','none');
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
%         set(gca,'YDir','reverse')
hold off
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
        savefig(gcf,[Directory 'Section ' OpenFiles '\Velocity_North.fig']);
    case 'East'
        title('Velocity East');
        [LowerVE,UpperVE] = limites(VelVector); 
        caxis([LowerVE UpperVE]);
        savefig(gcf,[Directory 'Section ' OpenFiles '/Velocity_East.fig']);
    case 'Upstream'
        title('Velocity Upstream');
        [LowerVU,UpperVU] = limites(VelVector); 
        caxis([LowerVU UpperVU]);
        savefig(gcf,[Directory 'Section ' OpenFiles '/Velocity_Upstream.fig']);
    case 'Error'
        title('Error Velocity');
        [LowerVD,UpperVD] = limites(VelVector); 
        caxis([LowerVD UpperVD]);
        savefig(gcf,[Directory 'Section ' OpenFiles '/Velocity_Upstream.fig']);
    case 'Magnitude'
        title('Velocity magnitude');
        [~,UpperVMag,StdMag] = limites(VelVector); 
        UpperVMag1 = ceil(max(max(VelVector))*10)/10;
        hold on
        savefig(gcf,[Directory 'Section ' OpenFiles '/Velocity_Magnitude.fig']);
    case 'Transverse'
        [LowerVR,UpperVR] = limites(VelVector); 
        caxis([LowerVR UpperVR]);
        title('Velocity transverse')
        savefig(gcf,[Directory 'Section ' OpenFiles '/Velocity_Transversal.fig']);
    case 'Longitudinal'
        title('Velocity Longitudinal')
        [LowerVL,UpperVL] = limites(VelVector); 
        caxis([LowerVL UpperVL]);
%         Ylim = get(gca,'Ylim');
%         subplot(1,10,10)
%         [VelLPerfAve, VelLPerfAveDepth]= VelPerfil(A.System.Cell_Size,A.System.Cell_Start,Sumnum,VelL,NumOfCellsExit,xx);
%         [VelUPerfAve, VelUPerfAveDepth]= VelPerfil(A.System.Cell_Size,A.System.Cell_Start,Sumnum,VelU,NumOfCellsExit,xx);
%         VelUPerfAveDepthx = zeros(size(VelLPerfAveDepth));
%         quiverc(VelUPerfAveDepthx,VelLPerfAveDepth,VelLPerfAve,VelUPerfAve)
%         Xlim2 = get(gca,'Xlim');
%         set(gca,'YDir','reverse')
%         ylim([0,-min(depth)+.1*diff([0,-min(depth)])]);
%         set(gca,'ytick',[]);set(gca,'xtick',[])
%         text(Xlim2(2),mean(Ylim),'Velocidade Up','horizontalAlignment','center','verticalAlignment','bottom','Rotation',-90)
        savefig(gcf,[Directory 'Section ' OpenFiles '/Velocity_Longitudinal.fig']);
end

end