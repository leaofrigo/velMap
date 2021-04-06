function SideVortS(handles)

data = get(handles.Options2D.VortSide,'UserData');
NUMBOX = data(1)+1;
perc = data(2);
yygrid = handles.SideVort.yygrid; 
xxgrid = handles.SideVort.xxgrid; 
Vel = handles.SideVort.Vel;
track = handles.SideVort.track;
depth = handles.SideVort.depth;
StartEdge = 0;
global qh qq qpatch
disttot = max(track)-min(track);
distperc = disttot*perc;
h=figure('WindowStyle','docked');
ScreenSize = get(0,'ScreenSize');
ScreenSizeWidth = ScreenSize(3);
ScreenSizeHeight = ScreenSize(4);
slider = uicontrol('Style', 'slider',...
        'Min',min(track),'Max',max(track)-distperc,'Value',min(track),...
        'Position', round([130/1920*ScreenSizeWidth 650*ScreenSizeHeight/1080,...
        300/1920*ScreenSizeWidth 20*ScreenSizeHeight/1080]),...
        'Callback', @surfzlim); 
ax(1) = axes('Position',[0.25,0.1,0.7,0.7]);
ax(2) = axes('Position',[0.1,0.75,0.23,0.23]);
set(ax(2),'Ytick',[])
set(ax(2),'Xtick',[])
d=[];
columname = {'Variable','U_prime(m/s)','L(m)','Moment (m^2/s)'};
columnformat = {'char','numeric','numeric','numeric'};
table = uitable(h,'Data',d,'ColumnName',columname,'ColumnEditable',[false false false false],...
            'RowName',[],'ColumnFormat', columnformat,'ColumnWidth',{65/1920*ScreenSizeWidth,...
            75/1920*ScreenSizeWidth,70/1920*ScreenSizeWidth,85/1920*ScreenSizeWidth});
set(table,'Position',round([130/1920*ScreenSizeWidth 350*ScreenSizeHeight/1080,...
    300/1920*ScreenSizeWidth 300*ScreenSizeHeight/1080]))
        

p = pcolor(xxgrid,yygrid,Vel);
set(p,'EdgeColor','none')
[LowerVN,UpperVN] = limites(Vel);
caxis([LowerVN,UpperVN])
C = hsv(250);
colormap(flipud(C([1:4:16,14:50,50:2:end-125,end-125:end-50,end-50:4:end-30],:))) % make the color bar look like the River Servoyor Live
cc=colorbar('south');
set(get(cc,'xlabel'),'String', 'Velocidade(m/s)');
initpos = get(cc,'Position');
initfontsize = get(cc,'FontSize');
set(cc,'Position',[initpos(1)+initpos(3)/2 initpos(2) initpos(3)/2 initpos(4)/2],...
    'FontSize',initfontsize*.75)
hold on
plot(track,depth)
xlim([min(track),max(track)])
if StartEdge == 1
    set(ax(2),'XDir','reverse')
end

quiverplot(min(track));
set(gcf,'toolbar','figure');
axes(ax(2))
p1 = patch('Faces',[1 2 3 4],'Vertices',[min(track) min(depth);...
    min(track) 0;min(track)+distperc 0;min(track)+distperc min(depth)],'FaceColor','none');
    
    function surfzlim(source,callbackdata)
        val = get(source,'Value');
        delete(p1)
        axes(ax(2))
        p1 = patch('Faces',[1 2 3 4],'Vertices',[val min(depth);...
            val 0;val+distperc 0;val+distperc min(depth)],'FaceColor','none');
%         delete(qh)
%         delete(qq)
%         delete(qpatch)
        quiverplot(val);


        
    end

    function quiverplot(val)
        

        TFx = and(xxgrid>val,xxgrid<val+distperc);      
        yTF1 = yygrid(TFx);
        Veltemp = Vel(TFx);
        TFNotNan = ~isnan(Veltemp);
       
        miny = min(yTF1(TFNotNan));
        maxy = max(yTF1(TFNotNan));      
        boxes = linspace(maxy,miny,NUMBOX);
        boxwidth = distperc;
        boxhight = abs(boxes(2)-boxes(1));
%         boxarea = boxwidth*boxhight;
        
        TFbox = zeros(length(yTF1),length(boxes)-1);
        VelAve= zeros(1,length(boxes)-1);
        yloc =VelAve;
        
%         uprime=yloc;
        for j = 1:length(boxes)-1
            if j ==length(boxes)-1
                TFbox(:,j) = and(yTF1(:,1)<=boxes(j),yTF1(:,1)>=boxes(j+1)); 
            else
                TFbox(:,j) = and(yTF1(:,1)<=boxes(j),yTF1(:,1)>boxes(j+1)); 
            end
            VelAve(j) = nanmean(Veltemp(logical(TFbox(:,j))));
            yloc(j) = nanmean(boxes(j:j+1));
%             uprime1(j) = Qave(j)/boxarea;
        end
        uprime1 = VelAve;
        ubar=nanmean(VelAve);
        uprime = zeros(1,length(VelAve)-1);
        un = cell(1,length(VelAve)-1);
        L=uprime;
        for kk = 1:length(uprime)
            uprime(kk) = VelAve(kk)-VelAve(kk+1);
            L(kk) = yloc(kk)-yloc(kk+1);
            un{kk} = ['U_' num2str(kk)];
        end
        Mom = uprime.*L; %uprime = u_n-u_n+1 L=dist(u_n,u_n+1) ubar=mean(u_1...u_n)
        dd1 = [{'U_bar'} num2cell(ubar) {''} {''}];

        dd = [un' num2cell(uprime)' num2cell(L)' num2cell(Mom)'];
        dd= [dd1;dd];
        set(table,'Data',dd)
        axes(ax(1))
        qh = quiver(ax(1),ones(1,length(boxes)-1)*nanmean([val val+distperc]),...
            yloc,VelAve,zeros(1,length(boxes)-1),'b','AutoScale','Off');
        xlimi=get(ax(1),'Xlim');
        ylimi=get(ax(1),'Ylim');
        ylabel('Depth(m)')
        xlabel('Location(m)')
        hold on
        scal = nanmean(uprime1)/5;
        scal = abs(scal); %in case of negative numbers
        n=0;
        while (scal*10^n<1)
            n=n+1;
        end
        scal =round(scal*10^n)/10^n;
%         qq = quiver(ax(1),xlimi(2)-(xlimi(2)-xlimi(1))*.025,ylimi(1)+(ylimi(2)-ylimi(1))*.1,1,0,0.1,'b');
        qq = quiver(ax(1),xlimi(2)-(xlimi(2)-xlimi(1))*.025,ylimi(1)+(ylimi(2)-ylimi(1))*.1,scal,...
            0,'b','AutoScale','Off');
        text(xlimi(2)-(xlimi(2)-xlimi(1))*.0225,ylimi(1)+(ylimi(2)-ylimi(1))*.12,[num2str(scal) ' m/s']);
%         scalefactor = get(qq,'AutoScaleFactor');
        qpatch = patch('Faces',[1 2 3 4],'Vertices',[xlimi(2)-(xlimi(2)-xlimi(1))*.025...
            ylimi(1)+(ylimi(2)-ylimi(1))*.05;...
            xlimi(2)-(xlimi(2)-xlimi(1))*.025 ylimi(1)+(ylimi(2)-ylimi(1))*.15;...
            xlimi(2)-(xlimi(2)-xlimi(1))*.025+scal ylimi(1)+(ylimi(2)-ylimi(1))*.15;...
            xlimi(2)-(xlimi(2)-xlimi(1))*.025+scal ylimi(1)+(ylimi(2)-ylimi(1))*.05],'FaceColor','none');
        hold off
        axes(ax(2))
        set(cc,'Position',[initpos(1)+initpos(3)/2 initpos(2) initpos(3)/2 initpos(4)/2],...
            'FontSize',initfontsize*.75)
    end

end
    
    