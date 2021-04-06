function sideVort(handles)

data = get(handles.Options2D.VortSide,'UserData');
NUMBOX = data(1)+1;
perc = data(2);
Face = handles.SideVort.Face;
VerticeXXYY=handles.SideVort.VerticeXXYY;
Vel = handles.SideVort.Vel;
track = handles.SideVort.track;
depth = handles.SideVort.depth;
StartEdge = handles.SideVort.StartEdge;
ft=fittype('a*(exp(x*b))');
global qh qq qpatch
disttot = max(track)-min(track);
distperc = disttot*perc;
h2 = figure('WindowStyle','normal');
h=figure('WindowStyle','normal');
set(h,'units','normalized','outerposition',[0 0 1 1])
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

p = patch('Faces',Face,'Vertices',VerticeXXYY);
set(p,'FaceColor','flat',...
'FaceVertexCData',Vel,...
'CDataMapping','scaled','EdgeColor','none');
[LowerVN,UpperVN] = limites(Vel);
caxis([LowerVN,UpperVN])
C = hsv(250);
colormap(flipud(C([1:4:16,14:50,50:2:end-125,end-125:end-50,end-50:4:end-30],:))) % make the color bar look like the River Servoyor Live
hold on
plot(track,depth)
xlim([min(track),max(track)])
if StartEdge == 1
    set(ax(2),'XDir','reverse')
end
x = Face(:,1:2);%h.x=x;
xx = VerticeXXYY(x(:,1),1);%h.xx=xx;
xx2 = VerticeXXYY(x(:,2),1);%h.xx2=xx2;
xxx=[xx xx2];%h.xxx=xxx;
y=Face(:,2:3);%h.y=y;
yy = VerticeXXYY(y(:,1),2);%h.yy=yy;
yy2 = VerticeXXYY(y(:,2),2);%h.yy2=yy2;
yyy = [yy yy2];%h.yyy=yyy;
yyymean=mean(yyy')';%h.yyymean=yyymean;
Face1 = Face(:,1);
quiverplot(min(track));
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
        

        TFx = and(xx>val,xx2<val+distperc);%h.TFx=TFx;        
        yTF1 = yy(TFx);%h.yTF1=yTF1;
        yTF2 = yy2(TFx);%h.yTF2=yTF2;
        yyTF = [yTF1 yTF2];%h.yyTF=yyTF;
        xTF1 = xx(TFx);%h.xTF1=xTF1;
        xTF2 = xx2(TFx);%h.xTF2=xTF2;
        xxTF = [xTF1 xTF2];%h.xxTF=xxTF;
        xyTF = [xxTF yyTF];%h.xyTF=xyTF;
        miny = min(yyTF(:));%h.miny=miny;
        maxy = max(yyTF(:));%h.maxy=maxy;
        boxes = linspace(maxy,miny,NUMBOX);%h.boxes=boxes;
        boxwidth = distperc;
        boxhight = abs(boxes(2)-boxes(1));
        boxarea = boxwidth*boxhight;
        TFbox = zeros(length(yyy),length(boxes)-1);
        TFface = cell(1,length(boxes)-1);
        Area = TFface;
        TFVel=Area;
        Areasum = 0;
        Qcell = Area;
        Qave= zeros(1,length(boxes)-1);
        yloc =Qave;  
        uprime=yloc;
        uprime1=yloc;
        for j = 1:length(boxes)-1
            if j ==length(boxes)-1
                TFbox(:,j) = and(yyymean(:,1)<=boxes(j),yyymean(:,1)>=boxes(j+1)); 
            else
                TFbox(:,j) = and(yyymean(:,1)<=boxes(j),yyymean(:,1)>boxes(j+1)); 
            end
            TFbox(:,j) = and(TFbox(:,j),TFx);
            TFface{j} = Face1(logical(TFbox(:,j)));
            TFVel{j} = Vel(logical(TFbox(:,j)));
            AreaTemp = zeros(length(TFface{j}),1);
            for k=1:length(TFface{j});
                f=TFface{j}(k);
                AreaTemp(k) = abs(VerticeXXYY(f+1,1)-VerticeXXYY(f,1))*abs(VerticeXXYY(f+2,2)-VerticeXXYY(f,2));
                Areasum = Areasum + sum(AreaTemp(k));
            end
            Qcell{j} = TFVel{j}.*AreaTemp;
            Area{j} = AreaTemp;
            Qave(j) = nanmean(Qcell{j});
            yloc(j) = nanmean(boxes(j:j+1));
            uprime1(j) = Qave(j)/boxarea;
        end
        ubar=nanmean(uprime1);
        uprime = zeros(1,length(uprime1)-1);
        L=uprime;
        un = cell(1,length(uprime1)-1);
        for kk = 1:length(uprime)
            uprime(kk) = uprime1(kk)-uprime1(kk+1);
            L(kk) = yloc(kk)-yloc(kk+1);
            un{kk} = ['U_' num2str(kk)];
        end
        Mom = uprime.*L;
        dd1 = [{'U_bar'} num2cell(ubar) {''} {''}];
        dd = [un' num2cell(uprime)' num2cell(L)' num2cell(Mom)'];
        dd= [dd1;dd];
        set(table,'Data',dd)
        axes(ax(1))
        qh = quiver(ax(1),ones(1,length(boxes)-1)*mean([val val+distperc]),...
                yloc,uprime1,zeros(1,length(boxes)-1),'b','AutoScale','Off');
        
        Depthmean = nanmean(depth(and(track>val,track<val+distperc)));
        DepthminMinusMax = abs(min(depth(and(track>val,track<val+distperc)))-max(depth(and(track>val,track<val+distperc))));
        DepthSTD = std(depth(and(track>val,track<val+distperc)));
        trye = [yloc;uprime1];
        trye(:,length(trye)+1)=[min(depth(and(track>val,track<val+distperc)));0];
        trye(1,:)=trye(1,:)-min(trye(1,:));
%         assignin('base', 'trye', trye)
%         assignin('base', 'Depthmean', Depthmean)
%         assignin('base', 'DepthminMinusMax', DepthminMinusMax)
%         assignin('base', 'DepthSTD', DepthSTD)
        f = fit(trye(2,:)',trye(1,:)',ft);
        coeff = coeffvalues(f);
        u_tau = 0.41/coeff(2);
        mu_turb = u_tau*1000*.41*linspace(0,max(trye(1,:)),10);
        figure(h2)
        clf
        fplot(f,[0,1.025*max(trye(2,:))]);
        hold on
        scatter(trye(2,:),trye(1,:));
        ylabel('Depth');xlabel('Velocity Magnitude Profile')
        
        figure(h)
        
        xlimi=get(ax(1),'Xlim');
        ylimi=get(ax(1),'Ylim');
        ylabel('Depth(m)')
        xlabel('Location(m)')
        hold on
        scal = mean(uprime1)/5;
        scal = abs(scal); %in case of negative numbers
        n=0;
        while (scal*10^n<1)
            n=n+1;
        end
        scal =round(scal*10^n)/10^n;
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
        
        
    end

end
    
    