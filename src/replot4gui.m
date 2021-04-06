function replot = replot4gui(handles)
Directory = handles.Directory;
OpenFiles = handles.OpenFiles;
Long = cell(length(OpenFiles),1);
Lat = Long;
Quiver3X = Lat;
Quiver3Y = Lat;
Quiver3U = Lat;
Quiver3V = Lat;
StartEdge = Lat;Quiver3VectorU1 =Lat;Quiver3VectorV1 = Lat; Quiver3VectorW1 = Lat;
Quiver3VectorX1 = Lat; Quiver3VectorY1 = Lat; Quiver3VectorZ1 = Lat; SumnumN = Lat;
for n=1:length(OpenFiles)
    A = load ([Directory OpenFiles{n}]);
    velocidades = A.WaterTrack.Velocity;
    Dist = A.Summary.Track;
    StartEdge{n} = A.Setup.startEdge;
    y = diff([0 0; Dist]); % calcular diferenca entre pontos
    clear('velE','velN','velU','velD');
    velE(:,:) = velocidades(:,1,:);
    velN(:,:) = velocidades(:,2,:);
    velU(:,:) = velocidades(:,3,:);
    velD(:,:) = velocidades(:,4,:);
    TamanhoMatrix = size(velocidades);
    xx = zeros(1,TamanhoMatrix(3)+1);
    for i = 1:TamanhoMatrix(3)
        xx(i+1) = xx(i) + norm(y(i,:)); 
    end
    %Calculando um linha reta entre primeiro e ultimo ponto
    TotalDist = xx(end) - xx(1);
    xxNorm = xx/TotalDist;      
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
    %inicializar/reiniciar variaveis
    AverageVelocityVector = zeros(TamanhoMatrix(3),3);
    zz = zeros(1,TamanhoMatrix(3)-1);
    AverageVelocityMag = zeros(1,TamanhoMatrix(3));
    Quiver3VectorX =[];Quiver3VectorY=[];Quiver3VectorZ =[];
    Long{n} = A.GPS.Longitude;
    Lat{n} = A.GPS.Latitude;
    Lattemp = A.GPS.Latitude;
    Longtemp = A.GPS.Longitude;
    Sumnum = [0 zz];
    Sumnum2=[0 zz];
    X = [];
    Y = [];
    cellDepths =[];
    cellDepthsCentro = [];
    for k = 1:TamanhoMatrix(3)
        [Sumnum(k+1),Sumnum2(k+1),X,Y,AverageVelocityVector(k,:),AverageVelocityMag(k),...
            cellDepths,zz(k+1),cellDepthsCentro,...
            ~,~,~,~,~,...
            Quiver3VectorX,Quiver3VectorY,Quiver3VectorZ,~] =...
            CriarVetoresPosVelAve(ComecoCell(k),CellSize(k),NumOfCells(k),cellDepths,...
            cellDepthsCentro,xx(k),velN(:,k),velE(:,k),velU(:,k),velD(:,k),Sumnum(k),Sumnum2(k),X,Y,...
            Quiver3VectorX,Quiver3VectorY,Quiver3VectorZ,Lattemp(k),Longtemp(k));        
    end
    Quiver3VectorU = velE(~isnan(velE));
    Quiver3VectorV = velN(~isnan(velN));
    Quiver3VectorW = velU(~isnan(velN));
    Sumnum = Sumnum+1;
    quiverX = zeros(length(Sumnum)-1,1);
    quiverY=quiverX;quiverU=quiverX;quiverV=quiverX;

    for w = 1:length(Sumnum)-1
        quiverX(w) = mean(Quiver3VectorX(Sumnum(w):Sumnum(w+1)-1));
        quiverY(w) = mean(Quiver3VectorY(Sumnum(w):Sumnum(w+1)-1));
        quiverU(w) = mean(Quiver3VectorU(Sumnum(w):Sumnum(w+1)-1));
        quiverV(w) = mean(Quiver3VectorV(Sumnum(w):Sumnum(w+1)-1));
    end
    
    SumnumN{n} = Sumnum;
    Quiver3VectorU1{n} = Quiver3VectorU;
    Quiver3VectorV1{n} = Quiver3VectorV;
    Quiver3VectorW1{n} = Quiver3VectorW;
    Quiver3VectorX1{n} = Quiver3VectorX;
    Quiver3VectorY1{n} = Quiver3VectorY;
    Quiver3VectorZ1{n} = Quiver3VectorZ; 
    
    
    Quiver3X{n} = quiverX;
    Quiver3Y{n} = quiverY;
    Quiver3U{n} = quiverU;
    Quiver3V{n} = quiverV;
end


replot.Lat = Lat;
replot.Long = Long;
replot.Quiver3X = Quiver3X;
replot.Quiver3Y = Quiver3Y;
replot.Quiver3U = Quiver3U;
replot.Quiver3V = Quiver3V;
replot.StartEdge = StartEdge;

replot.Quiver3allX = Quiver3VectorX1;
replot.Quiver3allY = Quiver3VectorY1;
replot.Quiver3allZ = Quiver3VectorZ1;
replot.Quiver3allU = Quiver3VectorU1;
replot.Quiver3allV = Quiver3VectorV1;
replot.Quiver3allW = Quiver3VectorW1;
replot.Sumnum = SumnumN;