function [] = simulationRiver3DEDDY(Quiver3VectorXTot,Quiver3VectorYTot,...
    Quiver3VectorZTot,EDDYVECTOR,...
    YrR,XrR,FrR,handles)
Space = get(handles.Options3D.NumofPointsBat,'UserData');
data = get(handles.Options3D.SimulOpt,'UserData');
pausetime = data(1);
VectorRotation = data(2:end);
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

VEL3D = griddata(Quiver3VectorXTot,Quiver3VectorYTot,Quiver3VectorZTot,EDDYVECTOR,X1,Y1,Z1);
% assignin('base', 'EDDYVECTOR', EDDYVECTOR)
% assignin('base', 'Quiver3VectorXTot', Quiver3VectorXTot)
% assignin('base', 'Quiver3VectorYTot', Quiver3VectorYTot)
% assignin('base', 'Quiver3VectorZTot', Quiver3VectorZTot)
% assignin('base', 'X1', X1)
% assignin('base', 'Y1', Y1)
% assignin('base', 'Z1', Z1)
% assignin('base', 'YrR', YrR)
% assignin('base', 'XrR', XrR)
% assignin('base', 'FrR', FrR)

h1 = figure('DeleteFcn',@func);
c=uicontrol(h1,'Style','pushbutton','String','Pause','Callback',@pausefunc);
cpos = get(c,'Position');
c1pos = cpos;
c1pos(1) = c1pos(1)+c1pos(3);
uicontrol(h1,'Style','pushbutton','String','Stop','Position',c1pos,'Callback',@stopfunc);
c2pos= c1pos;
c2pos(1) = c2pos(1)+c2pos(3);
uicontrol(h1,'Style','pushbutton','String','Run','Position',c2pos,'Callback',@runfunc);
c3pos = c2pos;
c3pos(1) = c1pos(1)-c1pos(3);
c3pos(3) = c1pos(3)*2;
c3pos(2) = c3pos(2)+c3pos(4);
slider = uicontrol(h1,'Style','slider','Min',1,'Max',Space,'Value',1,'Position',c3pos,'Callback',@slidefunc);
c4pos = c3pos;
c4pos(1) = c3pos(1)+c3pos(3);
c4pos(3) = c1pos(3);
uicontrol(h1,'Style','pushbutton','String','Pin','Position',c4pos,'Callback',@savesesfunc);
c5pos=cpos;
c5pos(2)=cpos(2)+2*cpos(4);
uicontrol(h1,'Style','pushbutton','String','Clear','Position',c5pos,'Callback',@clearfunc);

C=colormap;
Gr = surf(YrR,XrR,FrR);
colormap('gray')
freezeColors 
colormap(C)
set(Gr,'EdgeColor', 'none');
Xlim = get(gca,'Xlim'); Ylim =get(gca,'Ylim');Zlim = get(gca,'Zlim');
i=1;
j=1;
[~, UpperLim,s] = limites(VEL3D);
UpperLim = UpperLim-s;
ylabel('Latitude (deg)')
xlabel('Longitude (deg)')
zlabel('Depth (m)')

h2 = figure('Visible','Off');
global fe fe2 d
fe2=0;
d=1;
Vel = EliminateVSmallerThanBat(X1,Y1,Z1,VEL3D,XrR,YrR,FrR);
% Vel=VEL3D;
xd= zeros(Space,Space,Space);
yd=xd;zd=xd;
for k = 1:Space
    figure(h2)
    set(h2,'Visible','off')
    hslice = surf(ones(size(x1))*x1(k),y1,meshgrid(z1));
    rotate(hslice,VectorRotation(1:3),VectorRotation(end),[mean([xmin,xmax]) mean([ymin ymax])...
        mean([zmin zmax])])
    xd(:,:,k) = get(hslice,'XData');
    yd(:,:,k) = get(hslice,'YData');
    zd(:,:,k) = get(hslice,'ZData'); 
    
end
delete(h2)
% save('C:\Users\Roberta\Desktop\Ricardo\Variables Save\simul','xd','yd','zd','X1','Y1','Z1','Vel','VectorRotation')
play=true;
set( h1, 'toolbar', 'figure' )
playanimation
    function playanimation
        while play;

            figure(h1)

            hold on
            fe1 = slice(X1,Y1,Z1,Vel,xd(:,:,i),yd(:,:,i),zd(:,:,i));
        %     fe = slice(X1,Y1,Z1,Vel,0,Y1(i,1,1),0);
            if j == 1
                colorbar('eastoutside')
                j = 0;
                h_bar = findobj(gcf,'Tag','Colorbar');
                set(get(h_bar,'xlabel'),'String', 'Eddy Viscosity {\nu}_t (m^2.s^{-1})');
            end
            caxis([0, UpperLim]);
            set(fe1,'EdgeColor', 'none');
            axis([Xlim Ylim Zlim]);
            i=i+1;
            if i>Space
                i=1;
            end    
            pause(pausetime);
            delete(fe1)
            set(slider,'Value',i)


        end
    end
    function func(source,event)
        play=false;
    end
    function pausefunc(source,event)
        play=false;
        fe = slice(X1,Y1,Z1,Vel,xd(:,:,i),yd(:,:,i),zd(:,:,i));
        caxis([0, UpperLim]);
        set(fe,'EdgeColor', 'none');
        axis([Xlim Ylim Zlim]);
    end
    function stopfunc(source,event)        
        play=false;
        fe = slice(X1,Y1,Z1,Vel,xd(:,:,i),yd(:,:,i),zd(:,:,i));
        caxis([0, UpperLim]);
        set(fe,'EdgeColor', 'none');
        axis([Xlim Ylim Zlim]);
        i=1;
    end
    function runfunc(source,event)
        delete(fe)
        play=true;
        playanimation
    end
    function slidefunc(source,event)
        value = get(source,'Value');
        value = round(value);
        if value == 0 
            value =1;
        elseif value == Space+1;
            value =Space;
        end
        i=value;
        delete(fe)
        fe = slice(X1,Y1,Z1,Vel,xd(:,:,i),yd(:,:,i),zd(:,:,i));
        caxis([0, UpperLim]);
        set(fe,'EdgeColor', 'none');
        axis([Xlim Ylim Zlim]);
    end
    function savesesfunc(source,event)
        fe2(d) = slice(X1,Y1,Z1,Vel,xd(:,:,i),yd(:,:,i),zd(:,:,i));
        caxis([0, UpperLim]);
        set(fe2(d),'EdgeColor', 'none');
        axis([Xlim Ylim Zlim]);
        d=d+1;
    end
    function clearfunc(source,event)
        d=d-1;        
        for hh=1:d
            delete(fe2(hh));
        end
        d=1;
        fe2=0;
    end


end