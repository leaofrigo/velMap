clc; close all; clear all;
Directory = 'C:\Users\Roberta\Desktop\Ricardo\detalhamento do tauri adcp\Hold on\';
Bathymetry = 0; %Se Tiver batimetria real colocar 1, 0 interpola entre seçoes pra fazer batimetria em 3d;
% progressbar(0)
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
    TotalDist = xx(end) - xx(1);
    xxNorm = xx/TotalDist;      %Calculating a straight line between first and last point
    DistMag = sqrt((Dist(1,1)-Dist(end,1))^2+(Dist(1,2)-Dist(end,2))^2);
    xx = xxNorm*DistMag;
    
    ComecoCell = A.System.Cell_Start;
    CellSize = A.System.Cell_Size;
    NumOfCells  = A.Summary.Cells;
    AverageVelocityVector = zeros(TamanhoMatrix(3),3);
    zz = zeros(1,TamanhoMatrix(3)-1);
    AverageVelocityMag = zeros(1,TamanhoMatrix(3));
    Quiver3VectorX =[];Quiver3VectorY=[];Quiver3VectorZ =[];Quiver3VectorU=[];Quiver3VectorV =[];Quiver3VectorW=[];
    Long = A.GPS.Longitude;
    Lat = A.GPS.Latitude;
    Sumnum = [0 zz];
    Sumnum2=0;
    X = [];
    Y = [];
    cellDepths =[];
    for k = 1:TamanhoMatrix(3)
        temp = -ComecoCell(k):-CellSize(k):-NumOfCells(k)*CellSize(k)-ComecoCell(k);
        temp = temp';
        cellDepths = [cellDepths;temp];
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
        figure(1)
        %plot N
        hold on
        imagesc([((xx(k)+xx(k+1))/2-(xx(k+1)-xx(k))/6) ((xx(k)+xx(k+1))/2+(xx(k+1)-xx(k))/6)],...
            [0-(0+temp(1)*(1+1/(numcel-1))^-1/(2*(numcel-1))),-temp(1)+(0+temp(1)*(1+1/(numcel-1))^-1/(2*(numcel-1)))],TopVelNExt);   
        imagesc([((xx(k)+xx(k+1))/2-(xx(k+1)-xx(k))/6) ((xx(k)+xx(k+1))/2+(xx(k+1)-xx(k))/6)],...
            [-temp(end)-(-temp(end)+depth(k))*(1+1/(numcel-1))^-1/(2*(numcel-1)),...
            -depth(k)+(-temp(end)+depth(k))*(1+1/(numcel-1))^-1/(2*(numcel-1))],BottomVelNExt);
        imagesc([((xx(k)+xx(k+1))/2-(xx(k+1)-xx(k))/6) ((xx(k)+xx(k+1))/2+(xx(k+1)-xx(k))/6)],[ComecoCell(k),-temp(end)],VelNTemp);
        figure(2)
        %plot E
        hold on
        imagesc([((xx(k)+xx(k+1))/2-(xx(k+1)-xx(k))/6) ((xx(k)+xx(k+1))/2+(xx(k+1)-xx(k))/6)],...
            [0-(0+temp(1)*(1+1/(numcel-1))^-1/(2*(numcel-1))),-temp(1)+(0+temp(1)*(1+1/(numcel-1))^-1/(2*(numcel-1)))],TopVelEExt);   
        imagesc([((xx(k)+xx(k+1))/2-(xx(k+1)-xx(k))/6) ((xx(k)+xx(k+1))/2+(xx(k+1)-xx(k))/6)],...
            [-temp(end)-(-temp(end)+depth(k))*(1+1/(numcel-1))^-1/(2*(numcel-1)),...
            -depth(k)+(-temp(end)+depth(k))*(1+1/(numcel-1))^-1/(2*(numcel-1))],BottomVelEExt);
        imagesc([((xx(k)+xx(k+1))/2-(xx(k+1)-xx(k))/6) ((xx(k)+xx(k+1))/2+(xx(k+1)-xx(k))/6)],[ComecoCell(k),-temp(end)],VelETemp);
        figure(3)
        %plot U
        hold on
        imagesc([((xx(k)+xx(k+1))/2-(xx(k+1)-xx(k))/6) ((xx(k)+xx(k+1))/2+(xx(k+1)-xx(k))/6)],...
            [0-(0+temp(1)*(1+1/(numcel-1))^-1/(2*(numcel-1))),-temp(1)+(0+temp(1)*(1+1/(numcel-1))^-1/(2*(numcel-1)))],TopVelUExt);   
        imagesc([((xx(k)+xx(k+1))/2-(xx(k+1)-xx(k))/6) ((xx(k)+xx(k+1))/2+(xx(k+1)-xx(k))/6)],...
            [-temp(end)-(-temp(end)+depth(k))*(1+1/(numcel-1))^-1/(2*(numcel-1)),...
            -depth(k)+(-temp(end)+depth(k))*(1+1/(numcel-1))^-1/(2*(numcel-1))],BottomVelUExt);
        imagesc([((xx(k)+xx(k+1))/2-(xx(k+1)-xx(k))/6) ((xx(k)+xx(k+1))/2+(xx(k+1)-xx(k))/6)],[ComecoCell(k),-temp(end)],VelUTemp);
        figure(4)
        %plot Diff
        hold on
        imagesc([((xx(k)+xx(k+1))/2-(xx(k+1)-xx(k))/6) ((xx(k)+xx(k+1))/2+(xx(k+1)-xx(k))/6)],...
            [0-(0+temp(1)*(1+1/(numcel-1))^-1/(2*(numcel-1))),-temp(1)+(0+temp(1)*(1+1/(numcel-1))^-1/(2*(numcel-1)))],TopVelDExt);   
        imagesc([((xx(k)+xx(k+1))/2-(xx(k+1)-xx(k))/6) ((xx(k)+xx(k+1))/2+(xx(k+1)-xx(k))/6)],...
            [-temp(end)-(-temp(end)+depth(k))*(1+1/(numcel-1))^-1/(2*(numcel-1)),...
            -depth(k)+(-temp(end)+depth(k))*(1+1/(numcel-1))^-1/(2*(numcel-1))],BottomVelDExt);
        imagesc([((xx(k)+xx(k+1))/2-(xx(k+1)-xx(k))/6) ((xx(k)+xx(k+1))/2+(xx(k+1)-xx(k))/6)],[ComecoCell(k),-temp(end)],VelDTemp);
        figure(5)
        %plot Mag
        hold on
        imagesc([((xx(k)+xx(k+1))/2-(xx(k+1)-xx(k))/6) ((xx(k)+xx(k+1))/2+(xx(k+1)-xx(k))/6)],...
            [0-(0+temp(1)*(1+1/(numcel-1))^-1/(2*(numcel-1))),-temp(1)+(0+temp(1)*(1+1/(numcel-1))^-1/(2*(numcel-1)))],TopVelMagExt);   
        imagesc([((xx(k)+xx(k+1))/2-(xx(k+1)-xx(k))/6) ((xx(k)+xx(k+1))/2+(xx(k+1)-xx(k))/6)],...
            [-temp(end)-(-temp(end)+depth(k))*(1+1/(numcel-1))^-1/(2*(numcel-1)),...
            -depth(k)+(-temp(end)+depth(k))*(1+1/(numcel-1))^-1/(2*(numcel-1))],BottomVelMagExt);
        imagesc([((xx(k)+xx(k+1))/2-(xx(k+1)-xx(k))/6) ((xx(k)+xx(k+1))/2+(xx(k+1)-xx(k))/6)],[ComecoCell(k),-temp(end)],VelMagTemp);  
        if StartEdge ==1
            quiver((xx(k)+xx(k+1))/2,-temp(1:end-1),-VelETemp,-VelNTemp)

        else
            quiver((xx(k)+xx(k+1))/2,-temp(1:end-1),VelETemp,-VelNTemp)
        end        
        Quiver3VectorX = [Quiver3VectorX; ones(length(temp)-1,1)*A.GPS.Longitude(k)];
        Quiver3VectorY = [Quiver3VectorY; ones(length(temp)-1,1)*A.GPS.Latitude(k)];
        Quiver3VectorZ = [Quiver3VectorZ; temp(1:end-1)];
        Quiver3VectorU = [Quiver3VectorU; VelETemp];
        Quiver3VectorV = [Quiver3VectorV; VelNTemp];
        Quiver3VectorW = [Quiver3VectorW; VelUTemp];      
    end
    C = hsv(250);
    EdgeDist0 = A.Setup.Edges_0__DistanceToBank;EdgeDist1 = A.Setup.Edges_1__DistanceToBank;
    type0 = A.Setup.Edges_0__Method; type1 = A.Setup.Edges_1__Method;
    track = xx;
    [Area0,depthnew,tracknew] = edges(EdgeDist0,depth,type0,track,0);
    [Area1,depthnew1,tracknew1] = edges(EdgeDist1,depthnew,type1,tracknew,1);
    
    for kk=1:5;
        figure(kk)
        colorbar('south')
        h_bar = findobj(gcf,'Tag','Colorbar');
        set(get(h_bar,'xlabel'),'String', 'Velocidade(m/s)');
        initpos = get(h_bar,'Position');
        initfontsize = get(h_bar,'FontSize');
        set(h_bar,'Position',[initpos(1)+initpos(3)/2 initpos(2) initpos(3)/2 initpos(4)/2],...
            'FontSize',initfontsize*.75)
        plot(tracknew1, -depthnew1)
        colormap(flipud(C([1:4:16,14:50,50:2:end-125,end-125:end-50,end-50:4:end-30],:))) % make the color bar look like the River Servoyor Live
        xlim([tracknew1(1),tracknew1(end)])
        ylim([0,-min(depth)+.1*diff([0,-min(depth)])]);
        xlabel('Largura (m)')
        ylabel('Profundidade (m)')
        if StartEdge == 1
            set(gca,'XDir','reverse')
        end
        set(gca,'YDir','reverse')
        hold off
    end
    velMag = sqrt(velN.^2+velE.^2+velU.^2);
    figure(1)
    title('Velocidade na direção Norte');
    [LowerVN,UpperVN] = limites(velN); 
    caxis([LowerVN UpperVN]);
    set(gcf, 'Position', [680, 558, 560*3, 420*.75]);
    saveas(figure(1),[Directory 'Seção ' OpenFiles{n}(1:end-4) '\Velocidade_Norte.fig']);
    saveas(figure(1),[Directory 'Seção ' OpenFiles{n}(1:end-4) '\Velocidade_Norte.jpg']);
    figure(2)
    title('Velocidade na direção Leste');
    [LowerVE,UpperVE] = limites(velE); 
    caxis([LowerVE UpperVE]);
    set(gcf, 'Position', [680, 558, 560*3, 420*.75]);
    saveas(figure(2),[Directory 'Seção ' OpenFiles{n}(1:end-4) '/Velocidade_Leste.fig']);
    saveas(figure(2),[Directory 'Seção ' OpenFiles{n}(1:end-4) '/Velocidade_Leste.jpg']);
    figure(3)
    title('Velocidade para cima');
    [LowerVU,UpperVU] = limites(velU); 
    caxis([LowerVU UpperVU]);
    set(gcf, 'Position', [680, 558, 560*3, 420*.75]);
    saveas(figure(3),[Directory 'Seção ' OpenFiles{n}(1:end-4) '/Velocidade_Upstream.fig']);
    saveas(figure(3),[Directory 'Seção ' OpenFiles{n}(1:end-4) '/Velocidade_Upstream.jpg']);
    figure(4)
    title('Diferenca de velocidade (Erro de Velocidade)');
    [LowerVD,UpperVD] = limites(velD); 
    caxis([LowerVD UpperVD]);
    set(gcf, 'Position', [680, 558, 560*3, 420*.75]);
    saveas(figure(4),[Directory 'Seção ' OpenFiles{n}(1:end-4) '/Velocidade_Difference.fig']);
    saveas(figure(4),[Directory 'Seção ' OpenFiles{n}(1:end-4) '/Velocidade_Difference.jpg']);
    figure(5)
    title('Magnitude da velocidade');
    [LowerVMag,UpperVMag,StdMag] = limites(velMag); 
    UpperVMag1 = ceil(max(max(velMag))*10)/10;
    if UpperVMag1>floor(10*(UpperVMag+StdMag))/10
        caxis([0 floor(10*(UpperVMag+StdMag))/10]);
    else
        caxis([0 UpperVMag1])
    end
    set(gcf, 'Position', [680, 558, 560*3, 420*.75]);
    saveas(figure(5),[Directory 'Seção ' OpenFiles{n}(1:end-4) '/Velocidade_Magnitude.fig']);
    saveas(figure(5),[Directory 'Seção ' OpenFiles{n}(1:end-4) '/Velocidade_Magnitude.jpg']);
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
    VelL = zeros(1,Sumnum(end));
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
        figure(8)
        hold on
        imagesc([((xx(r)+xx(r+1))/2-(xx(r+1)-xx(r))/6) ((xx(r)+xx(r+1))/2+(xx(r+1)-xx(r))/6)],[ComecoCell(r),-temp1(end)],VelR(Sumnum(r)+1:Sumnum(r+1))');
        figure(9)
        hold on
        subplot(1,10,1:9)
        imagesc([((xx(r)+xx(r+1))/2-(xx(r+1)-xx(r))/6) ((xx(r)+xx(r+1))/2+(xx(r+1)-xx(r))/6)],[ComecoCell(r),-temp1(end)],VelL(Sumnum(r)+1:Sumnum(r+1))');
    end
%     [curlz,cav]= curl(X,Y,VelL,VelR);
    for rr=8:9;
        figure(rr)
        colorbar('south')
        h_bar = findobj(gcf,'Tag','Colorbar');
        set(get(h_bar,'xlabel'),'String', 'Velocidade(m/s)');
        initpos = get(h_bar,'Position');
        initfontsize = get(h_bar,'FontSize');
        set(h_bar,'Position',[initpos(1)+initpos(3)/2 initpos(2) initpos(3)/2 initpos(4)/2],...
            'FontSize',initfontsize*.75)
        plot(tracknew1, -depthnew1)
        colormap(flipud(C([1:4:16,14:50,50:2:end-125,end-125:end-50,end-50:4:end-30],:))) % make the color bar look like the River Servoyor Live
        xlim([tracknew1(1),tracknew1(end)])
%         ylim([0,6])
        xlabel('Largura (m)')
        ylabel('Profundidade (m)')
        if StartEdge == 1
            set(gca,'XDir','reverse')
        end
        set(gca,'YDir','reverse')
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
    saveas(figure(9),[Directory 'Seção ' OpenFiles{n}(1:end-4) '/Velocidade_direção_Rio.fig']);
    saveas(figure(9),[Directory 'Seção ' OpenFiles{n}(1:end-4) '/Velocidade_direção_Rio.jpg']);
    figure(8)
    [LowerVR,UpperVR] = limites(VelR); 
    caxis([LowerVR UpperVR]);
    title('Velocidade na direção transversal da Media do Rio')
    saveas(figure(8),[Directory 'Seção ' OpenFiles{n}(1:end-4) '/Velocidade_Transversal_Rio.fig']);
    saveas(figure(8),[Directory 'Seção ' OpenFiles{n}(1:end-4) '/Velocidade_Transversal_Rio.jpg']);
    close(figure(8),figure(9))
    ylim([0,-min(depth)+.1*diff([0,-min(depth)])]);
    
        
%     figure(6)
%     hold on
    LatTot = [LatTot; Lat];
    LongTot = [LongTot; Long];
    AverageVelocityVectorTot1 = [AverageVelocityVectorTot1;AverageVelocityVector(:,1) ];
    AverageVelocityVectorTot2 = [AverageVelocityVectorTot2;AverageVelocityVector(:,2) ];
    AverageVelocityVectorTot3 = [AverageVelocityVectorTot3;AverageVelocityVector(:,3) ];
%     quiverc(A.GPS.Longitude(1:end),A.GPS.Latitude(1:end),AverageVelocityVector(:,1),AverageVelocityVector(:,2));
    
    figure(7)
    hold on
    scale= .001;
    quiver3(Quiver3VectorX,Quiver3VectorY,Quiver3VectorZ,Quiver3VectorU,Quiver3VectorV,Quiver3VectorW,scale)
    close(figure(1),figure(2),figure(3),figure(4),figure(5))
    AreaBody = abs(FindCrossAreaRiver(depth,xx));
    AreaTot(n) = Area0+Area1+AreaBody;
    Perimeter(n) = FindPerimeterRiver(depthnew1,tracknew1);
    HidrRad(n) = AreaTot(n)/Perimeter(n);
    PercetageProcessed = n/length(OpenFiles)
end

QTotal = AreaTot.*TotalAverageVelocity';
htest = figure('units','normalized','outerposition',[0 0 1 1]);
h_bar = colorbar('North');
% colormap(h_bar,flipud(C([1:4:16,14:50,50:2:end-125,end-125:end-50,end-50:4:end-30],:)))
colormap(jet(64))
[LowerLimQuiver, UpperLimQuiver,sQuiver] = limites(sqrt(AverageVelocityVectorTot1.^2 + AverageVelocityVectorTot2.^2));
caxis([LowerLimQuiver, UpperLimQuiver]);
% quiverc(LongTot,LatTot,AverageVelocityVectorTot1,AverageVelocityVectorTot2);
hold on
quiverc(LongTot,LatTot,AverageVelocityVectorTot1,AverageVelocityVectorTot2);

xlabel('Longitude(deg)');ylabel('Latitude(deg)');
hp = plot_google_map;
% h_bar = colorbar('North');
set(get(h_bar,'xlabel'),'String', 'Velocidade(m/s)');
initpos = get(h_bar,'Position');
initfontsize = get(h_bar,'FontSize');
set(h_bar,'Position',[initpos(1)+initpos(3)/2 initpos(2) initpos(3)/2 initpos(4)/2],...
    'FontSize',initfontsize*.75)

% colormap(h_bar,flipud(C([1:4:16,14:50,50:2:end-125,end-125:end-50,end-50:4:end-30],:))) % make the color bar look like the River Servoyor Live
Xlim = get(gca,'Xlim');
Ylim = get(gca,'Ylim');
ax(2)=axes;
% Xlim = get(gca,'Xlim');
% Ylim = get(gca,'Ylim');

if Bathymetry == 0
    [Xr,Yr,Fr] = griddata(LongTot,LatTot,depth1,unique(LongTot),unique(LatTot)');
    [Xr1,Yr1,Fr1] = griddata(LongTot,LatTot,depth1,[LongTot(NumpFile(1:end-1)) LongTot(NumpFile(2:end)-1)],[LongTot(NumpFile(1:end-1)) LatTot(NumpFile(2:end)-1)]);
    IN = inpolygon(Xr,Yr,Xlim,Ylim);
    Fr(~IN) = NaN;
else
    [Xr,Yr,Fr] = batimetriaToc(min(LatTot),min(LongTot),max(LatTot),max(LongTot));
end
[C2,h2] = contourf(Xr,Yr,Fr);
alpha(0.1)
set(ax(2), 'XAxisLocation','bottom',...
             'YAxisLocation','right',...
             'Color','none');  
set(gca,'ytick',[]);set(gca,'xtick',[])
ylim(ax(2),Ylim)
xlim(ax(2),Xlim)
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
% linkaxes([h_bar,hc,hp],'xy')
% xlabel('Longitude(deg)');ylabel('Latitude(deg)');


saveas(figure(htest),[Directory 'Google Earth Files\Todas Seçoes.jpg']);
saveas(figure(htest),[Directory 'Google Earth Files\Todas Seçoes.fig']);
figure(7)
ylabel('Latitude (deg)')
xlabel('Longitude (deg)')
zlabel('Depth (m)')
surf(Xr,Yr,Fr);
colormap('gray')
saveas(figure(7),[Directory 'Google Earth Files\QuiverAll3D.fig']);
saveas(figure(7),[Directory 'Google Earth Files\QuiverAll3D.jpg']);
copyfile('C:\Users\Roberta\Desktop\Ricardo\googleearth\data\redcone.dae',[Directory 'Google Earth Files\redcone.dae'])
% progressbar(1)
% GoogleEarthQuiver([Directory 'Google Earth Files\Todas Seçoes'],LongTot,LatTot,AverageVelocityVectorTot1,...
%             AverageVelocityVectorTot2,AverageVelocityVectorTot3); 
close all



