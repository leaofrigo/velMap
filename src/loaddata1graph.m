function [Long, Lat,AverageVelocityVector,excel,handles] = loaddata1graph(Directory,Sec,...
    graf,smoo,Plot,quiv,handles)
A = load ([Directory Sec]);
velocidades = A.WaterTrack.Velocity;
depth = - A.Summary.Depth;
Dist = A.Summary.Track;
y = diff([0 0; Dist]); % calcular diferenca entre pontos
y1= sum(abs(y)); % somar a distancia absoluta entre pontos
largura=norm(y1); % Achar a magnitude
TamanhoMatrix = size(velocidades);
StartEdge = A.Setup.startEdge;
xx = zeros(1,TamanhoMatrix(3)+1);
n=1;
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
NumSmooth = get(handles.Options2D.Smoo2dOption,'UserData');
NumOfCellsY = NumSmooth(1)*max(NumOfCells);
NumOfCellsX = NumSmooth(2)*TamanhoMatrix(3);
%Calcular centro da celula em relacao a X
qq=0;
for ii = 1:length(xx)-1
    xxCentro(qq+1:qq+NumOfCells(ii)) = ones(1,NumOfCells(ii))*(mean([xx(ii) xx(ii+1)]));
    qq = NumOfCells(ii) + qq;
end
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
velE(:,:) = velocidades(:,1,:);
velN(:,:) = velocidades(:,2,:);
velU(:,:) = velocidades(:,3,:);
velD(:,:) = velocidades(:,4,:);
velMag = sqrt(velE.^2+velN.^2+velU.^2);
%%Go Through each column of data
Z0Ustar = get(handles.Options2D.Z0Ustar,'UserData');
ExtrapStruct.NumCellsExtrap = get(handles.Options2D.NumCellsExtrap,'UserData');
ExtrapStruct.expsetting.Z0 = Z0Ustar(1);
ExtrapStruct.expsetting.Ustar = Z0Ustar(2);
ExtrapStruct.power = Z0Ustar(3);
for k = 1:TamanhoMatrix(3)
    [Sumnum(k+1),Sumnum2(k+1),X,Y,AverageVelocityVector(k,:),AverageVelocityMag(k),...
        cellDepths,zz(k+1),cellDepthsCentro,...
        VelETemp,VelNTemp,VelUTemp,VelDTemp,VelMagTemp,...
        Quiver3VectorX,Quiver3VectorY,Quiver3VectorZ,temp] =...
        CriarVetoresPosVelAve(ComecoCell(k),CellSize(k),NumOfCells(k),cellDepths,...
        cellDepthsCentro,xx(k),velN(:,k),velE(:,k),velU(:,k),velD(:,k),Sumnum(k),Sumnum2(k),X,Y,...
        Quiver3VectorX,Quiver3VectorY,Quiver3VectorZ,Lat(k),Long(k));
    ExtrapStruct.VelVec = [VelETemp,VelNTemp,VelUTemp,VelDTemp,VelMagTemp];        
    if isempty(VelMagTemp) == 1
        continue
    end
    ExtrapStruct.Depth= depth(k);        
    ExtrapStruct.CellStart = ComecoCell(k);        
    ExtrapStruct.NumCellsInColumn=NumOfCells(k);
    ExtrapStruct.CellSize= CellSize(k);

    ExtrapStruct.metodo = A.Setup.extrapolation_Top_nFitType;
    ExtrapStruct.type = 'Top';
    [VelocityTopExtrp(:,k),Ytop(:,k)] = ExtrapolationVelocity(ExtrapStruct);

    ExtrapStruct.metodo = A.Setup.extrapolation_Bottom_nFitType;
    ExtrapStruct.type = 'Bottom';
    [VelocityBotExtrp(:,k),Ybot(:,k)] = ExtrapolationVelocity(ExtrapStruct);
    
end
VelExtap.VelocityTopExtrp = VelocityTopExtrp;
VelExtap.VelocityBotExtrp =VelocityBotExtrp;
VelExtap.Ytop=Ytop;
VelExtap.Ybot=Ybot;

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
[Area0,depthnew,tracknew] = edges(EdgeDist0,depth,type0,track,0,handles);
[Area1,depthnew1,tracknew1] = edges(EdgeDist1,depthnew,type1,tracknew,1,handles);
[AVERAGES] = CalcularMedias(xx,zz,AverageVelocityVector,AverageVelocityMag,VelExtap);
if type0 == 2
    QArea0 = Area0*0.3535*EdgeDist0*(-depth(1))*nanmean(velMag(:,1));
elseif type0 ==1
    QArea0 = 0.911*Area0*EdgeDist0*(-depth(1))*nanmean(velMag(:,1));
else QArea0=0;
end
if type1 == 2
    QArea1 = Area1*0.3535*EdgeDist1*(-depth(end))*nanmean(velMag(:,end));
elseif type1 ==1
    QArea1 = 0.911*Area1*EdgeDist1*(-depth(end))*nanmean(velMag(:,end));
else QArea1=0;
end

AverageQBody=AVERAGES.AverageQBody;
TotalAverageVelocityVector=AVERAGES.TotalAverageVelocityVector;
TotalAverageVelocity2 = AVERAGES.TotalAverageVelocity2;
TotalAverageVelocity = AVERAGES.TotalAverageVelocity;
AverageBothAverages=AVERAGES.AverageBothAverages;
AreaMeasured = AVERAGES.AreaMeasured;
Qtop = AVERAGES.QTop;
Qbot = AVERAGES.QBot;
Qbody = AreaMeasured*AverageBothAverages;
QMiddle = AVERAGES.QMiddle;

c(n)= -TotalAverageVelocityVector(n,2)/TotalAverageVelocityVector(n,1);
th(n) = atan(c(n))+pi;
VelL = zeros(1,Sumnum(end));
VelR = VelL;
VelU = VelR;
for r=1:TamanhoMatrix(3)
    [VelU,VelL,VelR,temp1]= CalcularVelLVelR(ComecoCell(r),CellSize(r),NumOfCells(r),...
        velN(:,r),velE(:,r),velU(:,r),Sumnum(r:r+1),th(n),VelU,VelL,VelR);      
end
VelRL = [VelR',VelL'];
VelRL = VelRL(2:end,:);
switch graf
    case 'North'
        Vel = Quiver3VectorV;      
    case 'East'
        Vel = Quiver3VectorU;        
    case 'Upstream'
        Vel = Quiver3VectorW;        
    case 'Magnitude'
        Vel = Quiver3VectorMag;
    case 'Error'
        Vel = Quiver3VectorErr; 
    case 'Transverse'
        Vel = VelRL(:,1);
    case 'Longitudinal'     
        Vel = VelRL(:,2);
end
if Plot ==1
    figure
    if smoo == 1
        p = SmoothFull(xxCentro,NumOfCellsX,cellDepths,cellDepthsCentro,NumOfCellsY,Vel,...
            Sumnum2,xx,graf,C,StartEdge,tracknew1,depthnew1,Sec(1:end-4),...
            XXQuiver,YYQuiver,Quiver3VectorU,Quiver3VectorV,Quiver3VectorW,Quiver3VectorMag,...
            Quiver3VectorErr,VelRL,Directory,quiv,handles);
        
        
    else
        p = CriarESalvarFigura(Face,VerticeXXYY,Vel,tracknew1,...
            depthnew1,graf,Directory,Sec(1:end-4),XXQuiver,YYQuiver,Quiver3VectorU,...
            Quiver3VectorV,Quiver3VectorW,Quiver3VectorMag,Quiver3VectorErr,VelRL,C,...
            StartEdge,quiv,handles);
    end
end

handles.SideVort.Face =Face;
handles.SideVort.VerticeXXYY =VerticeXXYY;
handles.SideVort.Vel =Vel;
handles.SideVort.track =tracknew1;
handles.SideVort.depth =depthnew1;
handles.SideVort.StartEdge =StartEdge;

handles.Eddy.Face =Face;
handles.Eddy.VerticeXXYY =VerticeXXYY;
handles.Eddy.Vel =Vel;
handles.Eddy.track =track(2:end);
handles.Eddy.depth =depth;
handles.Eddy.StartEdge =StartEdge;
handles.Eddy.Lat =Lat;
handles.Eddy.Long =Long;

if handles.exp2grid2D.V==1
    zone = handles.exp2grid2D.zone;
    pathgrid = handles.exp2grid2D.path;
    filenamegrid = handles.exp2grid2D.filename;
    [Vr1,Vr2,Vr3,Lat1,Lon1,utmx, utmy,zone,TFTot] = LoadDataToGridSingle(Long,Lat,zone,...
        AverageVelocityVector(:,1),AverageVelocityVector(:,2),AverageVelocityVector(:,3),...
        filenamegrid, pathgrid);
    idx = ~isnan(Vr1) & ~isnan(Vr2);
    Longitudeidx = Lon1(idx);
    Latitudeidx = Lat1(idx);
    utmxidx = utmx(idx);
    utmyidx = utmy(idx);
    VelocityEidx = Vr1(idx);
    VelocityNidx = Vr2(idx);
    VelocityUidx = Vr3(idx);
    
    minX = ones(1,length(Long));minXI = minX;
    for i=1:length(Long)
        minXtemp = abs(Long(i) - Longitudeidx);
        minYtemp = abs(Lat(i) - Latitudeidx);
        [disttemp, IDX]= min(sqrt(minXtemp.^2 + minYtemp.^2));
        minX(i) = disttemp;
        minXI(i) = IDX;
    end
    GridVar2D.Longitudeidx = Longitudeidx(minXI);
    GridVar2D.Latitudeidx = Latitudeidx(minXI);
    GridVar2D.VelocityEidx = VelocityEidx(minXI);
    GridVar2D.VelocityNidx = VelocityNidx(minXI);
    GridVar2D.VelocityUidx = VelocityUidx(minXI);
    GridVar2D.utmxidx = utmxidx(minXI);
    GridVar2D.utmyidx = utmyidx(minXI);
    
    GridVar2D.Longitude = Lon1;
    GridVar2D.Latitude = Lat1;
    GridVar2D.utmx = utmx;
    GridVar2D.utmy = utmy;
    GridVar2D.zone = zone;
    GridVar2D.VelocityE = Vr1;
    GridVar2D.VelocityN = Vr2;
    GridVar2D.VelocityU = Vr3;    
    GridVar2D.FirstIdx = and(TFTot(:),idx(:));
    GridVar2D.SecondIdx = minXI;
    save([Directory 'Section ' Sec(1:end-4) '\GridDataSec' Sec(1:end-4) '.mat'],'GridVar2D')
end

    
    
    
if handles.plotpreview2d == 1 || handles.topVorticity == 1 
    
    quiverX = [];
    quiverY=quiverX;quiverU=quiverX;quiverV=quiverX;quiv=[];
    VelL=VelL';
    for w = 1:length(Sumnum)-1        
        quiv(w) = mean(VelL(Sumnum(w):Sumnum(w+1)-1));
        quiverX(w) = mean(Quiver3VectorX(Sumnum(w):Sumnum(w+1)-1));
        quiverY(w) = mean(Quiver3VectorY(Sumnum(w):Sumnum(w+1)-1));
        quiverU(w) = mean(Quiver3VectorU(Sumnum(w):Sumnum(w+1)-1));
        quiverV(w) = mean(Quiver3VectorV(Sumnum(w):Sumnum(w+1)-1));
    end

    if handles.plotpreview2d == 1
        plotprev = figure;
        hold on
        pl = plot(Long,Lat,'b'); 
        q(n)=quiver(quiverX,quiverY,quiverU,quiverV,'b');
        set(q(n),'ShowArrowHead','off')
        legend(Sec)
        xlabel('Longitude(deg)')
        ylabel('Latitude(deg)')
        title('Plot Preview')
    end
    if handles.topVorticity == 1 
        plotprevlong = figure;
        hold on
        pl = plot(Long,Lat,'b'); 
        quivU=-quiv*cos(th(n));
        quivV=-quiv*sin(th(n));
        q1(n)=quiver(quiverX,quiverY,quivU,quivV,'b');
        set(q1(n),'ShowArrowHead','off')
        legend(Sec)
        xlabel('Longitude(deg)')
        ylabel('Latitude(deg)')
        title('Plot Preview Longitugional')
        uprime = zeros(length(quiv)-1,1);
        L=uprime;
        for upr=1:length(uprime)
            uprime(upr) = VelL(upr)-VelL(upr+1);
            L(upr) = sqrt((Dist(upr,1)-Dist(upr+1,1))^2+(Dist(upr,2)-Dist(upr+1,2))^2);
%             L(upr) = sqrt((quiverX(upr)-quiverX(upr+1))^2+(quiverY(upr)-quiverY(upr)));
        end
        ubar= mean(quiv);
        Mom = uprime.*L;
        
    end

end


AreaBody = abs(FindCrossAreaRiver(depth,xx));
AreaTot(n) = Area0+Area1+AreaBody;
Perimeter(n) = FindPerimeterRiver(depthnew1,tracknew1);
HidrRad(n) = AreaTot(n)/Perimeter(n);
excel(1).section = Sec;
excel(1).Total_Q = Qbody+Qtop+Qbot+QArea0+QArea1; 
excel(1).Ave_speed = AverageBothAverages;
excel(1).totaldist = TotalDist;
excel(1).diffxx1xxend = diff([xx(1) xx(end)]);
excel(1).maxdepthstart = max(A.System.Cell_Start);
excel(1).perimeter = Perimeter;
excel(1).AreaTot = AreaTot;
excel(1).hidrrad = HidrRad;
excel(1).AreaMeas = AreaMeasured;













