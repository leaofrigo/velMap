clc; close all; clear all;tic
Directory = 'C:\Users\Roberta\Desktop\Ricardo\detalhamento do tauri adcp\seçoes pra rodar\Novas\';
Bathymetry = 1; %Se Tiver batimetria real colocar 1, 0 interpola entre seçoes pra fazer batimetria em 3d;
Grid =1; %Fazer figura comparando com modelo(Necessita de arquivo),
FazerIm3d =0;
listing = dir(Directory);
a = cell(length(listing),1);
tf = zeros(length(listing),1);
k = a;
NumOfCellsExit = 10;
for i = 1:length(listing)
    a{i} = listing(i,1).name;
    k{i} = strfind(a{i}, '.mat');
    tf(i) = ~isempty(k{i});

end
OpenFiles = a(logical(tf));
mkdir([Directory 'Google Earth Files'])
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
Quiver3VectorXTot =[];Quiver3VectorYTot=[];Quiver3VectorZTot=[];Quiver3VectorUTot=[];Quiver3VectorVTot=[];Quiver3VectorWTot=[];
SumNumTot =[];
for n=1:length(OpenFiles)
    tic
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


    vel_E = velocidades(:,1,:);
    vel_N = velocidades(:,2,:);
    vel_U = velocidades(:,3,:);
    vel_D = velocidades(:,4,:);

    TamanhoMatrix = size(velocidades);
    velE = zeros(TamanhoMatrix(1),TamanhoMatrix(3));
    velN=velE;
    velU=velE;
    velD=velE;
    StartEdge = A.Setup.startEdge;
    xx = zeros(1,TamanhoMatrix(3)+1);
    NumpFile(n+1) = NumpFile(n)+TamanhoMatrix(3);

    for i = 1:TamanhoMatrix(3)
        velE(:,i) = vel_E(:,1,i);
        velN(:,i) = vel_N(:,1,i);
        velU(:,i) = vel_U(:,1,i);
        velD(:,i) = vel_D(:,1,i);
        xx(i+1) = xx(i) + norm(y(i,:)); 
    end
    TFvel = ~isnan(velD);
    VelECurl = velE(TFvel);
    VelNCurl = velN(TFvel);
    VelUCurl = velU(TFvel);
    
    TotalDist = xx(end) - xx(1);
    xxNorm = xx/TotalDist;      %Calculating a straight line between first and last point
    DistMag = sqrt((Dist(1,1)-Dist(end,1))^2+(Dist(1,2)-Dist(end,2))^2);
    xx = xxNorm*DistMag;
    xxCentro = zeros(1,length(xx)-1);  
    ComecoCell = A.System.Cell_Start;
    CellSize = A.System.Cell_Size;
    NumOfCells  = A.Summary.Cells;
    qq=0;
    for ii = 1:length(xx)-1
        xxCentro(qq+1:qq+NumOfCells(ii)) = ones(1,NumOfCells(ii))*(mean([xx(ii) xx(ii+1)]));
        qq = NumOfCells(ii) + qq;
    end
    AverageVelocityVector = zeros(TamanhoMatrix(3),3);
    zz = zeros(1,TamanhoMatrix(3)-1);
    AverageVelocityMag = zeros(1,TamanhoMatrix(3));
    Quiver3VectorX =[];Quiver3VectorY=[];Quiver3VectorZ =[];Quiver3VectorU=[];Quiver3VectorV =[];Quiver3VectorW=[];Quiver3VectorErr=[];
    Quiver3VectorMag =[];
    Long = A.GPS.Longitude;
    Lat = A.GPS.Latitude;
    Sumnum = [0 zz];
    Sumnum2=0;
    X = [];
    Y = [];
    cellDepths =[];
    cellDepthsCentro = [];
    for k = 1:TamanhoMatrix(3)
        temp = -ComecoCell(k):-CellSize(k):-NumOfCells(k)*CellSize(k)-ComecoCell(k);
        temp = temp';
        cellDepths = [cellDepths;temp];
        tempCentro = zeros(1,length(temp)-1);
        for jj = 1:length(temp)-1;
            tempCentro(jj) = mean([temp(jj) temp(jj+1)]);
        end
        cellDepthsCentro = [cellDepthsCentro,tempCentro];
        zz(k) = abs(diff([temp(1),temp(end)]));
%         xx(k+1) = xx(k) + norm(y(k,:));
        AveLat = mean(Lat);
        AveLong = mean(Long);
        VelNTemp = velN(:,k);
        VelNTemp = VelNTemp(~isnan(VelNTemp));
        VelETemp = velE(:,k);
        VelETemp = VelETemp(~isnan(VelETemp));
        VelUTemp = velU(:,k);
        VelUTemp = VelUTemp(~isnan(VelUTemp));        
        VelDTemp = velD(:,k);
        VelDTemp = VelDTemp(~isnan(VelDTemp));
        VelMagTemp = sqrt(VelNTemp.^2+VelETemp.^2+VelUTemp.^2);
        Sumnum(k+1) = Sumnum(k)+length(VelDTemp);
        Sumnum2(k+1) = Sumnum2(k)+length(temp);
        X = [X xx(k)*ones(1,length(VelDTemp))];
        Y = [Y;temp(1:end-1)];
        AverageVelocityVector(k,:) = [mean(VelETemp) mean(VelNTemp) mean(VelUTemp)];
        AverageVelocityMag(k) = mean(VelMagTemp);
                 
        [TopVelMagExt,TopVelEExt,TopVelNExt,TopVelUExt,TopVelDExt,BottomVelMagExt,BottomVelEExt,BottomVelNExt,BottomVelUExt,BottomVelDExt,numcel]...
            = TopBottomExtrapolation(A.Setup.extrapolation_Top_nFitType,...
            A.Setup.extrapolation_Top_nEntireProfil,A.Setup.extrapolation_Top_dExponent,...
            A.Setup.extrapolation_Top_nCells,VelMagTemp,VelETemp,VelNTemp,VelUTemp,VelDTemp,temp,depth(k));
  
        Quiver3VectorX = [Quiver3VectorX; ones(length(tempCentro),1)*A.GPS.Longitude(k)];
        Quiver3VectorY = [Quiver3VectorY; ones(length(tempCentro),1)*A.GPS.Latitude(k)];
        Quiver3VectorZ = [Quiver3VectorZ; tempCentro(1:end)];
        Quiver3VectorU = [Quiver3VectorU; VelETemp];
        Quiver3VectorV = [Quiver3VectorV; VelNTemp];
        Quiver3VectorW = [Quiver3VectorW; VelUTemp]; 
        Quiver3VectorErr = [Quiver3VectorErr; VelDTemp]; 
        Quiver3VectorMag = [Quiver3VectorMag; VelMagTemp]; 
    end
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
    
    
    for kk=1:5;
        figure(kk)
        p(kk) = patch('Faces',Face,'Vertices',VerticeXXYY);
        colorbar('south')
        set(p(kk),'FaceColor','flat',...
        'FaceVertexCData',Quiver3VectorVel(:,kk),...
        'CDataMapping','scaled','EdgeColor','none');
        h_bar = findobj(gcf,'Tag','Colorbar');
        set(get(h_bar,'xlabel'),'String', 'Velocidade(m/s)');
        initpos = get(h_bar,'Position');
        initfontsize = get(h_bar,'FontSize');
        set(h_bar,'Position',[initpos(1)+initpos(3)/2 initpos(2) initpos(3)/2 initpos(4)/2],...
            'FontSize',initfontsize*.75)
        hold on
        set(gcf, 'Position', [680, 558, 560*3, 420*.75]);
        plot(tracknew1, depthnew1)
        colormap(flipud(C([1:4:16,14:50,50:2:end-125,end-125:end-50,end-50:4:end-30],:))) % make the color bar look like the River Servoyor Live
        xlim([tracknew1(1),tracknew1(end)])
        ylim([min(depth)+.1*diff([0,min(depth)]),0]);
        xlabel('Largura (m)')
        ylabel('Profundidade (m)')
        if StartEdge == 1
            set(gca,'XDir','reverse')
        end
%         set(gca,'YDir','reverse')
        hold off
    end
    velMag = sqrt(velN.^2+velE.^2+velU.^2);
    figure(1)
    title('Velocidade na direção Norte');
    [LowerVN,UpperVN] = limites(velN); 
    caxis([LowerVN UpperVN]);
    savefig(figure(1),[Directory 'Seção ' OpenFiles{n}(1:end-4) '\Velocidade_Norte.fig']);
    clf
%     saveas(figure(1),[Directory 'Seção ' OpenFiles{n}(1:end-4) '\Velocidade_Norte.jpg']);
    figure(2)
    title('Velocidade na direção Leste');
    [LowerVE,UpperVE] = limites(velE); 
    caxis([LowerVE UpperVE]);
    savefig(figure(2),[Directory 'Seção ' OpenFiles{n}(1:end-4) '/Velocidade_Leste.fig']);
%     saveas(figure(2),[Directory 'Seção ' OpenFiles{n}(1:end-4) '/Velocidade_Leste.jpg']);
    clf
    figure(3)
    title('Velocidade para cima');
    [LowerVU,UpperVU] = limites(velU); 
    caxis([LowerVU UpperVU]);
    savefig(figure(3),[Directory 'Seção ' OpenFiles{n}(1:end-4) '/Velocidade_Upstream.fig']);
%     saveas(figure(3),[Directory 'Seção ' OpenFiles{n}(1:end-4) '/Velocidade_Upstream.jpg']);
    clf
    figure(4)
    title('Diferenca de velocidade (Erro de Velocidade)');
    [LowerVD,UpperVD] = limites(velD); 
    caxis([LowerVD UpperVD]);
    savefig(figure(4),[Directory 'Seção ' OpenFiles{n}(1:end-4) '/Velocidade_Difference.fig']);
%     saveas(figure(4),[Directory 'Seção ' OpenFiles{n}(1:end-4) '/Velocidade_Difference.jpg']);
    clf
    figure(5)
    title('Magnitude da velocidade');
    [LowerVMag,UpperVMag,StdMag] = limites(velMag); 
    UpperVMag1 = ceil(max(max(velMag))*10)/10;
    hold on
    if StartEdge ==1
        quiver(XXQuiver,YYQuiver,-Quiver3VectorU,Quiver3VectorV,.1)
    else
        quiver(XXQuiver,YYQuiver,Quiver3VectorU,Quiver3VectorV,.1)
    end
    if UpperVMag1>floor(10*(UpperVMag+StdMag))/10
        caxis([0 floor(10*(UpperVMag+StdMag))/10]);
    else
        caxis([0 UpperVMag1])
    end
    savefig(figure(5),[Directory 'Seção ' OpenFiles{n}(1:end-4) '/Velocidade_Magnitude.fig']);
    clf
%     saveas(figure(5),[Directory 'Seção ' OpenFiles{n}(1:end-4) '/Velocidade_Magnitude.jpg']);
%     AveLat = mean(Lat);
%     AveLong = mean(Long);
%     if length(Long) == length(AverageVelocityVector)
%        GoogleEarthQuiver(FileName,Long,Lat,AverageVelocityVector(:,1),...
%             AverageVelocityVector(:,2),AverageVelocityVector(:,3)); 
%     else
%     
%         GoogleEarthQuiver(FileName,Long(1:end-1),Lat(1:end-1),AverageVelocityVector(:,1),...
%             AverageVelocityVector(:,2),AverageVelocityVector(:,3));
%     end
    
    
    AverageQ(n,:) = (sum([xx(2:end)'.*zz'.*AverageVelocityVector(:,1) xx(2:end)'.*zz'.*AverageVelocityVector(:,2) ...
    xx(2:end)'.*zz'.*AverageVelocityVector(:,3)]));
    TotalAverageVelocityVector(n,:) = AverageQ(n,:)/sum(sum(zz'.*xx(2:end)'));
    TotalAverageVelocity2(n) = sum(xx(2:end).*zz.*AverageVelocityMag)/sum(sum(zz'.*xx(2:end)'));
    TotalAverageVelocity(n) = sqrt(sum(TotalAverageVelocityVector(n,:).^2));
    AverageBothAverages(n) = mean([TotalAverageVelocity(n),TotalAverageVelocity2(n)]);
    c(n)= -TotalAverageVelocityVector(n,2)/TotalAverageVelocityVector(n,1);
    th(n) = atan(c(n));
%     RotMat(:,:,n) = [cos(th(n)) - sin(th(n));sin(th(n)) cos(th(n))];
    VelL = zeros(1,Sumnum(end)-1);
    VelR = VelL;
    VelU = VelR;
    
    for r=1:TamanhoMatrix(3)-1
        temp1 = -ComecoCell(r):-CellSize(r):-NumOfCells(r)*CellSize(r)-ComecoCell(r);
        temp1 = temp1';
        VelNTemp = velN(:,r);
        VelNTemp = VelNTemp(~isnan(VelNTemp));
        VelETemp = velE(:,r);
        VelETemp = VelETemp(~isnan(VelETemp));  
        VelUTemp = velU(:,r);
        VelUTemp = VelUTemp(~isnan(VelUTemp));
        VelU(Sumnum(r)+1:Sumnum(r+1)) = VelUTemp;
        VelL(Sumnum(r)+1:Sumnum(r+1)) = VelETemp*cos(th(n)) - VelNTemp*sin(th(n));
        VelR(Sumnum(r)+1:Sumnum(r+1)) = VelETemp*sin(th(n)) + VelNTemp*cos(th(n));          
%         figure(8)
%         hold on
%         imagesc([((xx(r)+xx(r+1))/2-(xx(r+1)-xx(r))/6) ((xx(r)+xx(r+1))/2+(xx(r+1)-xx(r))/6)],[ComecoCell(r),-temp1(end)],VelR(Sumnum(r)+1:Sumnum(r+1))');
%         figure(9)
%         hold on
%         subplot(1,10,1:9)
%         imagesc([((xx(r)+xx(r+1))/2-(xx(r+1)-xx(r))/6) ((xx(r)+xx(r+1))/2+(xx(r+1)-xx(r))/6)],[ComecoCell(r),-temp1(end)],VelL(Sumnum(r)+1:Sumnum(r+1))');
    end
    VelRL = [VelR',VelL'];
    VelRL = VelRL(2:end,:);
    pr=zeros(2,1);
%     [curlz,cav]= curl(X,Y,VelL,VelR);
    for rr=8:9;
        figure(rr)
        if rr==9
            subplot(1,10,1:9)
        end
        pr(rr-7) = patch('Faces',Face,'Vertices',VerticeXXYY);
        colorbar('south')
        set(pr(rr-7),'FaceColor','flat',...
        'FaceVertexCData',VelRL(:,rr-7),...
        'CDataMapping','scaled','EdgeColor','none');
        h_bar = findobj(gcf,'Tag','Colorbar');
        set(get(h_bar,'xlabel'),'String', 'Velocidade(m/s)');
        initpos = get(h_bar,'Position');
        initfontsize = get(h_bar,'FontSize');
        set(h_bar,'Position',[initpos(1)+initpos(3)/2 initpos(2) initpos(3)/2 initpos(4)/2],...
            'FontSize',initfontsize*.75)
        hold on
        plot(tracknew1, depthnew1)
        colormap(flipud(C([1:4:16,14:50,50:2:end-125,end-125:end-50,end-50:4:end-30],:))) % make the color bar look like the River Servoyor Live
        xlim([tracknew1(1),tracknew1(end)])
        ylim([min(depth)+.1*diff([0,min(depth)]),0]);
        xlabel('Largura (m)')
        ylabel('Profundidade (m)')
        if StartEdge == 1
            set(gca,'XDir','reverse')
        end
%         set(gca,'YDir','reverse')
        hold off
    end
    title('Velocidade na direção da Media do Rio')
    [LowerVL,UpperVL] = limites(VelL); 
    caxis([LowerVL UpperVL]);
    Ylim = get(gca,'Ylim');
    subplot(1,10,10)
    [VelLPerfAve, VelLPerfAveDepth]= VelPerfil(A.System.Cell_Size,A.System.Cell_Start,Sumnum,VelL,NumOfCellsExit,xx);
    [VelUPerfAve, VelUPerfAveDepth]= VelPerfil(A.System.Cell_Size,A.System.Cell_Start,Sumnum,VelU,NumOfCellsExit,xx);
    VelUPerfAveDepthx = zeros(size(VelLPerfAveDepth));
    quiverc(VelUPerfAveDepthx,VelLPerfAveDepth,VelLPerfAve,VelUPerfAve)
    Xlim2 = get(gca,'Xlim');
    set(gca,'YDir','reverse')
    ylim([0,-min(depth)+.1*diff([0,-min(depth)])]);
    set(gca,'ytick',[]);set(gca,'xtick',[])
    text(Xlim2(2),mean(Ylim),'Velocidade Up','horizontalAlignment','center','verticalAlignment','bottom','Rotation',-90)
    savefig(figure(9),[Directory 'Seção ' OpenFiles{n}(1:end-4) '/Velocidade_direção_Rio.fig']);
%     saveas(figure(9),[Directory 'Seção ' OpenFiles{n}(1:end-4) '/Velocidade_direção_Rio.jpg']);
    figure(8)
    [LowerVR,UpperVR] = limites(VelR); 
    caxis([LowerVR UpperVR]);
    title('Velocidade na direção transversal da Media do Rio')
    savefig(figure(8),[Directory 'Seção ' OpenFiles{n}(1:end-4) '/Velocidade_Transversal_Rio.fig']);
%     saveas(figure(8),[Directory 'Seção ' OpenFiles{n}(1:end-4) '/Velocidade_Transversal_Rio.jpg']);
    close(figure(8),figure(9))
    ylim([0,-min(depth)+.1*diff([0,-min(depth)])]);

    LatTot = [LatTot; Lat];
    LongTot = [LongTot; Long];
    AverageVelocityVectorTot1 = [AverageVelocityVectorTot1;AverageVelocityVector(:,1) ];
    AverageVelocityVectorTot2 = [AverageVelocityVectorTot2;AverageVelocityVector(:,2) ];
    AverageVelocityVectorTot3 = [AverageVelocityVectorTot3;AverageVelocityVector(:,3) ];
    
%     figure(7)
%     hold on
%     scale= .001;
%     quiver3(Quiver3VectorX,Quiver3VectorY,Quiver3VectorZ,Quiver3VectorU,Quiver3VectorV,Quiver3VectorW,scale)
%     close(figure(1),figure(2),figure(3),figure(4),figure(5))
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
%     GoogleEarthQuiver([Directory 'Google Earth Files\Seção' OpenFiles{n}(1:end-4)],Long,Lat,AverageVelocityVector(:,1),...
%             AverageVelocityVector(:,2),AverageVelocityVector(:,3)); 
    PercetageProcessed = n/length(OpenFiles)
    toc
end


QTotal = AreaTot.*TotalAverageVelocity';
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
    % colormap(h_bar,flipud(C([1:4:16,14:50,50:2:end-125,end-125:end-50,end-50:4:end-30],:))) % make the color bar look like the River Servoyor Live
    % Xlim = get(gca,'Xlim');
    % Ylim = get(gca,'Ylim');
    % ax(2)=axes;
    % Xlim = get(gca,'Xlim');
    % Ylim = get(gca,'Ylim');

    if Bathymetry == 0
        [Xr,Yr,Fr] = griddata(LongTot,LatTot,depth1,unique(LongTot),unique(LatTot)');
    %     [Xr1,Yr1,Fr1] = griddata(LongTot,LatTot,depth1,[LongTot(NumpFile(1:end-1)) LongTot(NumpFile(2:end)-1)],[LongTot(NumpFile(1:end-1)) LatTot(NumpFile(2:end)-1)]);
        IN = inpolygon(Xr,Yr,[min(LongTot) max(LongTot)],[min(LatTot) max(LatTot)]);
        Fr(~IN) = NaN;
        set(ax(1), 'XAxisLocation','bottom',...
                 'YAxisLocation','left',...
                 'Color','none'); 
        [C2,h2] = contourf(Xr,Yr,Fr);
        alpha(0.1)

    %     set(gca,'ytick',[]);set(gca,'xtick',[])
    %     ylim(ax(2),Ylim)
    %     xlim(ax(2),Xlim)
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
    %     contourf(YrR,XrR,FrR-A.GPS.Altitude(1),'LineStyle','none');
        hold on
        set(ax(1), 'XAxisLocation','bottom',...
                     'YAxisLocation','left','Color','none');


    %     ylim(ax(2),Ylim)
    %     xlim(ax(2),Xlim)
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
    % saveas(htest2,[Directory 'Google Earth Files\QuiverAll3D.jpg']);
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
    % saveas(figure(htest),[Directory 'Google Earth Files\Todas Seçoes.jpg']);
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



