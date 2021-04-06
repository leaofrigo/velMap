clc; close all; clear all;
tic
%%Inputs
Bathymetry = 1; %Se Tiver batimetria real colocar 1,0 interpola entre seçoes pra fazer batimetria em 3d;
Grid =1; %Fazer figura comparando com modelo(Necessita de arquivo),
FazerIm3d =1;
NumOfCellsExit = 10; %Numero de flexas demostrando Vel upstream em Velocidade diraçao do rio
SmoothYN = 0 ; %0 se quiser demonstrar grafico como adcp leu, 1 se quizer interpolar valoes
Directory = 'C:\Users\Roberta\Desktop\Ricardo\detalhamento do tauri adcp\seçoes pra rodar\Novas\Sem Smooth\';
%% Pegar Arquivos para abrir / criar pasta
OpenFiles = DirectoryOpenFiles(Directory);
mkdir([Directory 'Google Earth Files'])
%% Inicializar Variaveis
AverageQ = zeros(length(OpenFiles),3);
TotalAverageVelocityVector = zeros(length(OpenFiles),3);
TotalAverageVelocity2 = zeros(length(OpenFiles),1);
TotalAverageVelocity = TotalAverageVelocity2;
AverageBothAverages = TotalAverageVelocity2;
depth1 =[];
LatTot = [];
LongTot=[];
AverageVelocityVectorTot1 =[];
AverageVelocityVectorTot2 = [];
AverageVelocityVectorTot3 = [];
c = zeros(1,length(OpenFiles));
th = c;
AreaTot = c;Perimeter = c;HidrRad =c;
NumpFile = ones(1,length(OpenFiles)+1);
Quiver3VectorXTot =[];Quiver3VectorYTot=[];Quiver3VectorZTot=[];Quiver3VectorUTot=[];
Quiver3VectorVTot=[];Quiver3VectorWTot=[];
SumNumTot =[];
%% Ir por Todos os arquivos
for n=1:length(OpenFiles)
     A = load ([Directory OpenFiles{n}]);
    velocidades = A.WaterTrack.Velocity;
    depth = - A.Summary.Depth;
    depth1 = [depth1;depth];
    mkdir([Directory 'Seção ' OpenFiles{n}(1:end-4)])    
    FileName = [Directory 'Google Earth Files\Seçao_' OpenFiles{n}(1:end-4)];
    %Calcular distancia viajada pelo barco
    Dist = A.Summary.Track;
    y = diff([0 0; Dist]); % calcular diferenca entre pontos
    y1= sum(abs(y)); % somar a distancia absoluta entre pontos
    largura=norm(y1); % Achar a magnitude


    velE(:,:) = velocidades(:,1,:);
    velN(:,:) = velocidades(:,2,:);
    velU(:,:) = velocidades(:,3,:);
    velD(:,:) = velocidades(:,4,:);

    TamanhoMatrix = size(velocidades);
    StartEdge = A.Setup.startEdge;
    xx = zeros(1,TamanhoMatrix(3)+1);
    NumpFile(n+1) = NumpFile(n)+TamanhoMatrix(3);

    for i = 1:TamanhoMatrix(3)
        xx(i+1) = xx(i) + norm(y(i,:)); 
    end
    TFvel = ~isnan(velD);
    VelECurl = velE(TFvel);
    VelNCurl = velN(TFvel);
    VelUCurl = velU(TFvel);
    %Calculando um linha reta entre primeiro e ultimo ponto
    TotalDist = xx(end) - xx(1);
    xxNorm = xx/TotalDist;      
    DistMag = sqrt((Dist(1,1)-Dist(end,1))^2+(Dist(1,2)-Dist(end,2))^2);
    xx = xxNorm*DistMag;
    xxCentro = zeros(1,length(xx)-1);  
    ComecoCell = A.System.Cell_Start;
    CellSize = A.System.Cell_Size;
    NumOfCells  = A.Summary.Cells;
    NumOfCellsY = 3*max(NumOfCells);
    NumOfCellsX = 3*TamanhoMatrix(3);
    %Calcular centro da celula em relacao a X
    qq=0;
    for ii = 1:length(xx)-1
        xxCentro(qq+1:qq+NumOfCells(ii)) = ones(1,NumOfCells(ii))*(mean([xx(ii) xx(ii+1)]));
        qq = NumOfCells(ii) + qq;
    end
    %inicializar/reiniciar variaveis
    AverageVelocityVector = zeros(TamanhoMatrix(3),3);
    zz = zeros(1,TamanhoMatrix(3)-1);
    AverageVelocityMag = zeros(1,TamanhoMatrix(3));
    Quiver3VectorX =[];Quiver3VectorY=[];Quiver3VectorZ =[];
    Long = A.GPS.Longitude;
    Lat = A.GPS.Latitude;
    Sumnum = [0 zz];
    Sumnum2=[0 zz];
    X = [];
    Y = [];
    cellDepths =[];
    cellDepthsCentro = [];
    AveLat = mean(Lat);
    AveLong = mean(Long);
    
    %%Go Through each column of data
    for k = 1:TamanhoMatrix(3)
        [Sumnum(k+1),Sumnum2(k+1),X,Y,AverageVelocityVector(k,:),AverageVelocityMag(k),...
            cellDepths,zz(k+1),cellDepthsCentro,...
            VelETemp,VelNTemp,VelUTemp,VelDTemp,VelMagTemp,...
            Quiver3VectorX,Quiver3VectorY,Quiver3VectorZ,temp] =...
            CriarVetoresPosVelAve(ComecoCell(k),CellSize(k),NumOfCells(k),cellDepths,...
            cellDepthsCentro,xx(k),velN(:,k),velE(:,k),velU(:,k),velD(:,k),Sumnum(k),Sumnum2(k),X,Y,...
            Quiver3VectorX,Quiver3VectorY,Quiver3VectorZ,Lat(k),Long(k));
        
        [TopVelMagExt,TopVelEExt,TopVelNExt,TopVelUExt,TopVelDExt,BottomVelMagExt,...
            BottomVelEExt,BottomVelNExt,BottomVelUExt,BottomVelDExt,numcel]...
            = TopBottomExtrapolation(A.Setup.extrapolation_Top_nFitType,...
            A.Setup.extrapolation_Top_nEntireProfil,A.Setup.extrapolation_Top_dExponent,...
            A.Setup.extrapolation_Top_nCells,VelMagTemp,VelETemp,VelNTemp,VelUTemp,VelDTemp,...
            temp,depth(k));        
        
    end
    Quiver3VectorW = velU(~isnan(velU));
    Quiver3VectorV = velN(~isnan(velN));
    Quiver3VectorU = velE(~isnan(velE));
    Quiver3VectorErr = velD(~isnan(velD));
    Quiver3VectorMag = sqrt(Quiver3VectorU.^2 + Quiver3VectorV.^2 + Quiver3VectorW.^2);
    Sumnum2=Sumnum2+1;
    Sumnum=Sumnum+1;
    Quiver3VectorPos = [Quiver3VectorX,Quiver3VectorY,Quiver3VectorZ];
    Quiver3VectorVel = [Quiver3VectorV,Quiver3VectorU,Quiver3VectorW,Quiver3VectorErr,Quiver3VectorMag];
    C = hsv(250);
    EdgeDist0 = A.Setup.Edges_0__DistanceToBank;EdgeDist1 = A.Setup.Edges_1__DistanceToBank;
    type0 = A.Setup.Edges_0__Method; type1 = A.Setup.Edges_1__Method;
    track = xx;
    [Face, VerticeXXYY,XXQuiver,YYQuiver] = CriarImagem(cellDepths,track,Sumnum,Sumnum2);
    [Area0,depthnew,tracknew] = edges(EdgeDist0,depth,type0,track,0);
    [Area1,depthnew1,tracknew1] = edges(EdgeDist1,depthnew,type1,tracknew,1);
    p=zeros(5,1);
    Direcao = [{'Norte'},{'Leste'},{'Cima'},{'Erro'},{'Magnitude'}];
    for kk=1:5;

        if SmoothYN == 1
            p(kk) = SmoothFull(xxCentro,NumOfCellsX,cellDepths,cellDepthsCentro,NumOfCellsY,Quiver3VectorVel(:,kk),...
                Sumnum2,xx,Direcao{kk},C,StartEdge,tracknew1,depthnew1,OpenFiles,n,...
                XXQuiver,YYQuiver,Quiver3VectorU,Quiver3VectorV,Directory);
        else
            p(kk) = CriarESalvarFigura(Face,VerticeXXYY,Quiver3VectorVel(:,kk),tracknew1,...
            depthnew1,Direcao{kk},Directory,OpenFiles,n,XXQuiver,YYQuiver,Quiver3VectorU,...
            Quiver3VectorV,C,StartEdge);
        end
    end
    
    [AverageQ(n,:),TotalAverageVelocityVector(n,:),TotalAverageVelocity2(n),...
        TotalAverageVelocity(n),AverageBothAverages(n)] = CalcularMedias(xx,zz,AverageVelocityVector,...
        AverageVelocityMag);
    c(n)= -TotalAverageVelocityVector(n,2)/TotalAverageVelocityVector(n,1);
    th(n) = atan(c(n));
    VelL = zeros(1,Sumnum(end));
    VelR = VelL;
    VelU = VelR;
    for r=1:TamanhoMatrix(3)
        [VelU,VelL,VelR,temp1]= CalcularVelLVelR(ComecoCell(r),CellSize(r),NumOfCells(r),...
            velN(:,r),velE(:,r),velU(:,r),Sumnum(r:r+1),th(n),VelU,VelL,VelR);      
    end
    VelRL = [VelR',VelL'];
    VelRL = VelRL(2:end,:);
    pr=zeros(2,1);
    Direcao1 =[{'Radial'} {'Longitugional'}];
    for rr= 1:2
        if SmoothYN == 1
            pr(rr) = SmoothFull(xxCentro,NumOfCellsX,cellDepths,cellDepthsCentro,NumOfCellsY,VelRL(:,rr),...
                Sumnum2,xx,Direcao1{rr},C,StartEdge,tracknew1,depthnew1,OpenFiles,n,...
                XXQuiver,YYQuiver,Quiver3VectorU,Quiver3VectorV,Directory);
        else
            pr(rr) = CriarESalvarFigura(Face,VerticeXXYY,VelRL(:,rr),tracknew1,...
                depthnew1,Direcao1{rr},Directory,OpenFiles,n,XXQuiver,YYQuiver,Quiver3VectorU,...
                Quiver3VectorV,C,StartEdge);
        end        
  
    end
    LatTot = [LatTot; Lat];
    LongTot = [LongTot; Long];
    AverageVelocityVectorTot1 = [AverageVelocityVectorTot1;AverageVelocityVector(:,1) ];
    AverageVelocityVectorTot2 = [AverageVelocityVectorTot2;AverageVelocityVector(:,2) ];
    AverageVelocityVectorTot3 = [AverageVelocityVectorTot3;AverageVelocityVector(:,3) ];
    Quiver3VectorXTot = [Quiver3VectorXTot;Quiver3VectorX];
    Quiver3VectorYTot = [Quiver3VectorYTot;Quiver3VectorY];
    Quiver3VectorZTot = [Quiver3VectorZTot;Quiver3VectorZ];
    Quiver3VectorUTot = [Quiver3VectorUTot;Quiver3VectorU];
    Quiver3VectorVTot = [Quiver3VectorVTot;Quiver3VectorV];
    Quiver3VectorWTot = [Quiver3VectorWTot;Quiver3VectorW];
    SumNumTot = [SumNumTot Sumnum];
    AreaBody = abs(FindCrossAreaRiver(depth,xx));
    AreaTot(n) = Area0+Area1+AreaBody;
    Perimeter(n) = FindPerimeterRiver(depthnew1,tracknew1);
    HidrRad(n) = AreaTot(n)/Perimeter(n);    
    PercetageProcessed = n/length(OpenFiles)
    toc
end
QTotal = AreaTot.*TotalAverageVelocity'; 
close all
if FazerIm3d == 1
    htest = figure;
    clf
    ax(3) = axes('Position',[0.25,0.1,0.7,0.7]);
    plot_google_map
    set(ax(3),'ytick',[]);set(ax(3),'xtick',[])
    ax(1) = axes('Position',[0.25,0.1,0.7,0.7]);
    [LowerLimQuiver, UpperLimQuiver,sQuiver] = limites(sqrt(AverageVelocityVectorTot1.^2 + AverageVelocityVectorTot2.^2));
    hold on
    xlabel('Longitude(deg)');ylabel('Latitude(deg)');

    if Bathymetry == 0
        [Xr,Yr,Fr] = griddata(LongTot,LatTot,depth1,unique(LongTot),unique(LatTot)');
        IN = inpolygon(Xr,Yr,[min(LongTot) max(LongTot)],[min(LatTot) max(LatTot)]);
        Fr(~IN) = NaN;
        set(ax(1), 'XAxisLocation','bottom',...
                 'YAxisLocation','left',...
                 'Color','none'); 
        [C2,h2] = contourf(Xr,Yr,Fr);
        alpha(0.1)

        for q=1:length(h2)
        set(h2(q),'LineStyle','none');
        end
        hc=colorbar('South');
        [LowerLimDepth, UpperLimDepth,sDepth] = limites(depth1);
        caxis([LowerLimDepth, 0])
        initposc = get(hc,'Position');
        initfontsizec = get(hc,'FontSize');
        set(hc,'Position',[initposc(1) initposc(2) initposc(3)/2 initposc(4)/2],...
            'FontSize',initfontsizec*.75)
        set(get(hc,'xlabel'),'String', 'Profundidade(m)');
        hp = plot_google_map;

        figure(htest2)
        hold on
        scale= .001;
        quiver3(Quiver3VectorXTot,Quiver3VectorYTot,Quiver3VectorZTot,Quiver3VectorUTot,Quiver3VectorVTot,Quiver3VectorWTot,scale)

        ylabel('Latitude (deg)')
        xlabel('Longitude (deg)')
        zlabel('Depth (m)')
        handle = surf(Xr,Yr,Fr);
        set(handle,'EdgeColor','none')
        colormap('gray')
    else
        [Xr,Yr,Fr] = batimetriaToc(min(LatTot),min(LongTot),max(LatTot),max(LongTot));
        [XrR,YrR] = meshgrid(linspace(min(LatTot),max(LatTot),100),linspace(min(LongTot),max(LongTot),100));    
        FrR = griddata(Xr,Yr,Fr,XrR,YrR);
        scatter(Yr,Xr,[],Fr,'filled')
        hold on
        set(ax(1), 'XAxisLocation','bottom',...
                     'YAxisLocation','left','Color','none');
        hc=colorbar('South');
        initposc = get(hc,'Position');
        initfontsizec = get(hc,'FontSize');
        set(hc,'Position',[initposc(1) initposc(2) initposc(3)/2 initposc(4)/2],...
            'FontSize',initfontsizec*.75)
        set(get(hc,'xlabel'),'String', 'Batimetria Altitude (m)');
        htest2 = figure
        scale= .001;
        quiver3(Quiver3VectorXTot,Quiver3VectorYTot,Quiver3VectorZTot,Quiver3VectorUTot,Quiver3VectorVTot,Quiver3VectorWTot,scale,'Color','g')
        hold on
        ylabel('Latitude (deg)')
        xlabel('Longitude (deg)')
        zlabel('Depth (m)')
        handle = surf(YrR,XrR,FrR-A.GPS.Altitude(1));
        set(handle,'EdgeColor','none')

    end
    savefig(htest2,[Directory 'Google Earth Files\QuiverAll3D.fig']);
    clf
    close(htest2)
    figure(htest)
    hold on
    ax(2) = axes('Position',[0.25,0.1,0.7,0.7]);
    quiverc(LongTot,LatTot,AverageVelocityVectorTot1,AverageVelocityVectorTot2,2);
    set(ax(2),'Color','none')
    set(ax(2),'ytick',[]);set(ax(2),'xtick',[])
    h_bar = colorbar('North');
    colormap(jet(64))
    set(get(h_bar,'xlabel'),'String', 'Velocidade(m/s)');
    initpos = get(h_bar,'Position');
    initfontsize = get(h_bar,'FontSize');
    set(h_bar,'Position',[initpos(1)+initpos(3)/2 initpos(2) initpos(3)/2 initpos(4)/2],...
        'FontSize',initfontsize*.75)
    caxis([LowerLimQuiver, UpperLimQuiver]);
    linkaxes([ax(1),ax(2),ax(3)],'xy');
    ylim(ax(1),'auto')
    Ylim=get(ax(1),'Ylim');Xlim=get(ax(1),'Xlim');
    ax(4) = axes('Position',[0.1,0.75,0.23,0.23]);
    box on
    set(ax(4),'ytick',[]);set(ax(4),'xtick',[])
    ylim(ax(4),[-33.75 5.266666666666667]);xlim(ax(4),[-73.98333333333333 -53.38333333333333])
    hold on
    plot_google_map
    patch([Xlim(1) Xlim(1) Xlim(2) Xlim(2)],[Ylim(1) Ylim(2) Ylim(1) Ylim(2)],zeros(4,1),'FaceColor','none','Marker','o')
    savefig(figure(htest),[Directory 'Google Earth Files\Todas Seçoes.fig']);
    clf
    copyfile('C:\Users\Roberta\Desktop\Ricardo\googleearth\data\redcone.dae',[Directory 'Google Earth Files\redcone.dae'])

    if Grid == 1 
        zone = '22 M';
        h = figure;
        [Vr1,Vr2,Vr3,Lat1,Lon1] = LoadDataToGrid(LongTot,LatTot,zone,AverageVelocityVectorTot1,AverageVelocityVectorTot2,AverageVelocityVectorTot3);
        clf
        ax1(3) = axes('Position',[0.25,0.1,0.7,0.7]);
        plot_google_map
        set(ax1(3),'ytick',[]);set(ax1(3),'xtick',[])
        ax1(1) = axes('Position',[0.25,0.1,0.7,0.7]);
        set(ax1(1), 'XAxisLocation','bottom',...
                 'YAxisLocation','left','Color','none');

        hold on
        xlabel('Longitude(deg)');ylabel('Latitude(deg)');
        if Bathymetry ==1
            scatter(Yr,Xr,[],Fr,'filled');
        else
            [C2,h2] = contourf(Xr,Yr,Fr,6,'LineStyle','none');
            alpha(0.1)
        end

        % contourf(Lon1,Lat1,sqrt(Vr1.^2 + Vr2.^2),6,'LineColor','none')
        hold on

        hc1=colorbar('South');
        initposc1 = get(hc1,'Position');
        initfontsizec1 = get(hc1,'FontSize');
        set(hc1,'Position',[initposc1(1) initposc1(2) initposc1(3)/2 initposc1(4)/2],...
            'FontSize',initfontsizec1*.75)
         if Bathymetry ==0
                set(get(hc1,'xlabel'),'String', 'Profundidade(m)');
                caxis([LowerLimDepth, 0])
         else
                set(get(hc1,'xlabel'),'String', 'Batimetria Altitude (m)');
         end

        ax1(2) = axes('Position',[0.25,0.1,0.7,0.7]);
        idx = ~isnan(Vr1) & ~isnan(Vr2); 
        quiverc(Lon1(idx),Lat1(idx),Vr1(idx),Vr2(idx));
        set(ax1(2),'Color','none')
        set(ax1(2),'ytick',[]);set(ax1(2),'xtick',[])
        h_bar1 = colorbar('North');
        colormap(jet(64))
        set(get(h_bar1,'xlabel'),'String', 'Velocidade(m/s)');
        initpos2 = get(h_bar1,'Position');
        initfontsize2 = get(h_bar1,'FontSize');
        set(h_bar1,'Position',[initpos2(1)+initpos2(3)/2 initpos2(2) initpos2(3)/2 initpos2(4)/2],...
            'FontSize',initfontsize2*.75)

        caxis([LowerLimQuiver, UpperLimQuiver]);
        linkaxes([ax1(1),ax1(2),ax1(3)],'xy');
        ylim(ax1(1),'auto')
        Ylim=get(ax1(1),'Ylim');Xlim=get(ax1(1),'Xlim');
        ax1(4) = axes('Position',[0.1,0.75,0.23,0.23]);
        box on
        set(ax1(4),'ytick',[]);set(ax1(4),'xtick',[])
        ylim(ax1(4),[-33.75 5.266666666666667]);xlim(ax1(4),[-73.9 -53.38333333333333])
        hold on
        plot_google_map
        patch([Xlim(1) Xlim(1) Xlim(2) Xlim(2)],[Ylim(1) Ylim(2) Ylim(1) Ylim(2)],zeros(4,1),'FaceColor','none','Marker','o')
        savefig(figure(htest),[Directory 'Google Earth Files\Todas Seçoes Grid.fig']);
    end
    close all
    GoogleEarthQuiver([Directory 'Google Earth Files\Todas Seçoes'],LongTot,LatTot,AverageVelocityVectorTot1,...
                AverageVelocityVectorTot2,AverageVelocityVectorTot3); 
else
    close all
end
toc





    




