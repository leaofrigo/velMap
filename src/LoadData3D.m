function [Quiver3VectorXTot,Quiver3VectorYTot,Quiver3VectorZTot,Quiver3VectorUTot,...
    Quiver3VectorVTot,Quiver3VectorWTot,YrR,XrR,FrR,alt,excel,GridVar,VarNme,VectorEddy]...
    = LoadData3D(Directory,batimetriav,OpenFiles,grid,PlSav,super,quiv3,kml3,smoo,simul,quiv,...
    filenamegrid, pathgrid,filenamebat, pathbat,zone,StreamLines,handles)
tic
if exist([Directory 'Analysis 3D'],'dir')
else
    mkdir([Directory 'Analysis 3D'])
end
%% Inicializar Variaveis
VarNme =[];
GridVar = [];
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
AreaTot = c;Perimeter = c;HidrRad =c;AreaMeasured=c;Qtop=c;Qbot=c;Qbody=c;QArea0=c;QArea1=c;
NumpFile = ones(1,length(OpenFiles)+1);
Quiver3VectorXTot =[];Quiver3VectorYTot=[];Quiver3VectorZTot=[];Quiver3VectorUTot=[];
Quiver3VectorVTot=[];Quiver3VectorWTot=[];
SumNumTot =[];
if handles.plotpreview == 1
    plotprev = figure;
    colours = distinguishable_colors(length(OpenFiles));
    legendstr =[];
end
%% Ir por Todos os arquivos
LatEddy=[];
LongEddy=[];
MUEddy = [];
VectorEddy =[];
DepthSimulEddy =[];
LatSimulEddy =[];
LongSimulEddy =[];
for n=1:length(OpenFiles)
    A = load ([Directory OpenFiles{n}]);
    velocidades = A.WaterTrack.Velocity;
    %Standard Deviation (error) used for EddyViscosity
    VelSTD = A.WaterTrack.Vel_StdDev;
    depth = - A.Summary.Depth;
    depth1 = [depth1;depth];
    if PlSav==1
        if exist([Directory 'Section ' OpenFiles{n}],'dir')
        else
            mkdir([Directory 'Section ' OpenFiles{n}])   
        end         
    end
    %Calcular distancia viajada pelo barco
    Dist = A.Summary.Track;
    y = diff([0 0; Dist]); % calcular diferenca entre pontos
    y1= sum(abs(y)); % somar a distancia absoluta entre pontos
    largura=norm(y1); % Achar a magnitude
    clear('velE','velN','velU','velD','VelSTDE','VelSTDN','VelSTDUp','VelSTDD');
    velE(:,:) = velocidades(:,1,:);
    velN(:,:) = velocidades(:,2,:);
    velU(:,:) = velocidades(:,3,:);
    velD(:,:) = velocidades(:,4,:);
    VelSTDE(:,:) = VelSTD(:,1,:);
    VelSTDN(:,:) = VelSTD(:,2,:);
    VelSTDUp(:,:) = VelSTD(:,3,:);
    VelSTDD(:,:) = VelSTD(:,4,:);
    
    
    
    velMag = sqrt(velE.^2+velN.^2+velU.^2);
    TamanhoMatrix = size(velocidades);
    StartEdge = A.Setup.startEdge;
    xx = zeros(1,TamanhoMatrix(3)+1);
    NumpFile(n+1) = NumpFile(n)+TamanhoMatrix(3);

    for i = 1:TamanhoMatrix(3)
        xx(i+1) = xx(i) + norm(y(i,:)); 
    end
    TFvel = ~isnan(velD);

    %Calculando um linha reta entre primeiro e ultimo ponto
    TotalDist = xx(end) - xx(1);
    xxNorm = xx/TotalDist;      
    DistMag = sqrt((Dist(1,1)-Dist(end,1))^2+(Dist(1,2)-Dist(end,2))^2);
    xx = xxNorm*DistMag;
    xxCentro = zeros(1,length(xx)-1);  
    ComecoCell = A.System.Cell_Start;
    CellSize = A.System.Cell_Size;
    NumOfCells  = A.Summary.Cells;
    NumSmooth = get(handles.Options3D.Smoo3dOption,'UserData');
    NumOfCellsY = NumSmooth(2)*max(NumOfCells);
    NumOfCellsX = NumSmooth(1)*TamanhoMatrix(3);
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
    
    Z0Ustar = get(handles.Options3D.Z0Ustar,'UserData');
    ExtrapStruct.NumCellsExtrap = get(handles.Options3D.NumCellsExtrap,'UserData');
    ExtrapStruct.expsetting.Z0 = Z0Ustar(1);
    ExtrapStruct.expsetting.Ustar = Z0Ustar(2);
    ExtrapStruct.power = Z0Ustar(3);
    %%Go Through each column of data
    VelocityTopExtrp = zeros(ExtrapStruct.NumCellsExtrap,TamanhoMatrix(3));
    VelocityBotExtrp =VelocityTopExtrp; Ytop =VelocityTopExtrp; Ybot=VelocityTopExtrp;
    for k = 1:TamanhoMatrix(3)
        [Sumnum(k+1),Sumnum2(k+1),X,Y,AverageVelocityVector(k,:),AverageVelocityMag(k),...
            cellDepths,zz(k+1),cellDepthsCentro,...
            VelETemp,VelNTemp,VelUTemp,VelDTemp,VelMagTemp,...
            Quiver3VectorX,Quiver3VectorY,Quiver3VectorZ,temp] =...
            CriarVetoresPosVelAve(ComecoCell(k),CellSize(k),NumOfCells(k),cellDepths,...
            cellDepthsCentro,xx(k),velN(:,k),velE(:,k),velU(:,k),velD(:,k),Sumnum(k),Sumnum2(k),X,Y,...
            Quiver3VectorX,Quiver3VectorY,Quiver3VectorZ,Lat(k),Long(k));
        
        ExtrapStruct.VelVec = VelMagTemp;       
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
%         
%         [TopVelMagExt,TopVelEExt,TopVelNExt,TopVelUExt,TopVelDExt,BottomVelMagExt,...
%             BottomVelEExt,BottomVelNExt,BottomVelUExt,BottomVelDExt,numcel]...
%             = TopBottomExtrapolation(A.Setup.extrapolation_Top_nFitType,...
%             A.Setup.extrapolation_Top_nEntireProfil,A.Setup.extrapolation_Top_dExponent,...
%             A.Setup.extrapolation_Top_nCells,VelMagTemp,VelETemp,VelNTemp,VelUTemp,VelDTemp,...
%             temp,depth(k));        
        
    end
    VelExtap.VelocityTopExtrp = VelocityTopExtrp;
    VelExtap.VelocityBotExtrp =VelocityBotExtrp;
    VelExtap.Ytop=Ytop;
    VelExtap.Ybot=Ybot;
    Quiver3VectorW = velU(~isnan(velU));
    Quiver3VectorV = velN(~isnan(velN));
    Quiver3VectorU = velE(~isnan(velE));
    Quiver3VectorErr = velD(~isnan(velD));
    %% Std Dev. Velocities - Eddy Viscosity
    VelSTDE = VelSTDE(~isnan(VelSTDE));
    VelSTDN = VelSTDN(~isnan(VelSTDN));
    VelSTDUp = VelSTDUp(~isnan(VelSTDUp));
    VelSTDD = VelSTDD(~isnan(VelSTDD));
    %    
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
        QArea0(n) = Area0*0.3535*EdgeDist0*(-depth(1))*nanmean(velMag(:,1));
    elseif type0 ==1
        QArea0(n) = 0.911*Area0*EdgeDist0*(-depth(1))*nanmean(velMag(:,1));
    else QArea0(n)=0;
    end
    if type1 == 2
        QArea1(n) = Area1*0.3535*EdgeDist1*(-depth(end))*nanmean(velMag(:,end));
    elseif type1 ==1
        QArea1(n) = 0.911*Area1*EdgeDist1*(-depth(end))*nanmean(velMag(:,end));
    else QArea1(n)=0;
    end
    
    AverageQBody=AVERAGES.AverageQBody;
    TotalAverageVelocityVector(n,:)=AVERAGES.TotalAverageVelocityVector;
    TotalAverageVelocity2(n,:) = AVERAGES.TotalAverageVelocity2;
    TotalAverageVelocity(n,:) = AVERAGES.TotalAverageVelocity;
    AverageBothAverages(n,:)=AVERAGES.AverageBothAverages;
    AreaMeasured(n) = AVERAGES.AreaMeasured;
    Qtop(n) = AVERAGES.QTop;
    Qbot(n) = AVERAGES.QBot;
    Qbody(n) = AreaMeasured(n)*AverageBothAverages(n);
    
    
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
    if PlSav == 1
        p=zeros(5,1);        
        Direcao = [{'North'},{'East'},{'Upstream'},{'Error'},{'Magnitude'}];
        for kk=1:5;  

            if smoo == 1
                p(kk) = SmoothFull(xxCentro,NumOfCellsX,cellDepths,cellDepthsCentro,NumOfCellsY,Quiver3VectorVel(:,kk),...
                    Sumnum2,xx,Direcao{kk},C,StartEdge,tracknew1,depthnew1,[OpenFiles{n}],...
                    XXQuiver,YYQuiver,Quiver3VectorU,Quiver3VectorV,Quiver3VectorW,...
                    Quiver3VectorMag,Quiver3VectorErr,VelRL,Directory,quiv,handles);
                clf
            else
                p(kk) = CriarESalvarFigura(Face,VerticeXXYY,Quiver3VectorVel(:,kk),tracknew1,...
                depthnew1,Direcao{kk},Directory,[OpenFiles{n}],XXQuiver,YYQuiver,Quiver3VectorU,...
                Quiver3VectorV,Quiver3VectorW,Quiver3VectorMag,Quiver3VectorErr,VelRL,C,StartEdge,quiv,handles);
                clf
            end
        end
        pr=zeros(2,1);
        Direcao1 =[{'Transverse'} {'Longitudinal'}];
        for rr= 1:2
            if smoo == 1
                pr(rr) = SmoothFull(xxCentro,NumOfCellsX,cellDepths,cellDepthsCentro,NumOfCellsY,VelRL(:,rr),...
                    Sumnum2,xx,Direcao1{rr},C,StartEdge,tracknew1,depthnew1,[OpenFiles{n}],...
                    XXQuiver,YYQuiver,Quiver3VectorU,Quiver3VectorV,Quiver3VectorW,...
                    Quiver3VectorMag,Quiver3VectorErr,VelRL,Directory,quiv,handles);
                clf
            else
                pr(rr) = CriarESalvarFigura(Face,VerticeXXYY,VelRL(:,rr),tracknew1,...
                    depthnew1,Direcao1{rr},Directory,[OpenFiles{n}],XXQuiver,YYQuiver,Quiver3VectorU,...
                    Quiver3VectorV,Quiver3VectorW,Quiver3VectorMag,Quiver3VectorErr,VelRL,C,StartEdge,...
                    quiv,handles);
                clf
            end  
           

        end

    end
    if handles.plotpreview == 1
        figure(plotprev);
        hold on
        pl(n) = plot(Long,Lat,'Color',colours(n,:)); 
        quiverX = [];
        quiverY=quiverX;quiverU=quiverX;quiverV=quiverX;UTMx=[];UTMy=[];UTMZone=[];varname=[];
        for w = 1:length(Sumnum)-1
            quiverX(w) = mean(Quiver3VectorX(Sumnum(w):Sumnum(w+1)-1));
            quiverY(w) = mean(Quiver3VectorY(Sumnum(w):Sumnum(w+1)-1));
            quiverU(w) = mean(Quiver3VectorU(Sumnum(w):Sumnum(w+1)-1));
            quiverV(w) = mean(Quiver3VectorV(Sumnum(w):Sumnum(w+1)-1));
        end
        q(n)=quiver(quiverX,quiverY,quiverU,quiverV,'Color',colours(n,:));
        set(q(n),'ShowArrowHead','off')
        legendstr = [legendstr {OpenFiles{n}}];
        TFY = ~isnan(quiverY);
        TFX = ~isnan(quiverX);
        quiverY = quiverY(and(TFY,TFX));
        quiverX= quiverX(and(TFY,TFX));
        quiverU= quiverU(and(TFY,TFX));
        quiverV= quiverV(and(TFY,TFX));
        VarNme{n} = genvarname(['LongitudeLatitudeMagSection_' OpenFiles{n}]);
        varname.Latitude = quiverY;
        varname.Longitude = quiverX;
        varname.MagnitudeVelocidade = sqrt(quiverU.^2+quiverV.^2);
        varname.VelocidadeLeste = quiverU;
        varname.VelocidadeNorte = quiverV;
        [UTMx,UTMy,UTMZone] = deg2utm(quiverY,quiverX);
        varname.UTMx = UTMx';
        varname.UTMy = UTMy';
        varname.UTMZone = UTMZone';
        assignin('caller', VarNme{n}, varname);
        if isdir([Directory 'Section ' OpenFiles{n}])==0
            mkdir([Directory 'Section ' OpenFiles{n}])
        end
        
%         save([Directory 'Seção ' OpenFiles{n} '/LongitudeLatitudeMagSeçao'],VarNme)
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
    excel(n).section = OpenFiles(n);
    excel(n).Total_Q = Qbody(n)+Qtop(n)+Qbot(n)+QArea0(n)+QArea1(n);    
    excel(n).Ave_speed = AverageBothAverages(n);
    excel(n).totaldist = TotalDist;
    excel(n).diffxx1xxend = diff([xx(1) xx(end)]);
    excel(n).maxdepthstart = max(A.System.Cell_Start);
    excel(n).perimeter = Perimeter(n);
    excel(n).AreaTot = AreaTot(n);
    excel(n).hidrrad = HidrRad(n);
    excel(n).AreaMeas = AreaMeasured(n);
    toc
    if handles.Eddy3DV == 1 || handles.EddySimulV == 1 || handles.EddyExpV ==1
        handles.Eddy.Face = Face;
        handles.Eddy.VerticeXXYY = VerticeXXYY;
        handles.Eddy.Vel = Quiver3VectorVel(:,5);
%         Changes requested by G.Thomas 2/8/16
        handles.Eddy.VelE = Quiver3VectorU;
        handles.Eddy.VelN = Quiver3VectorV;
        handles.Eddy.VelUp = Quiver3VectorW;
        handles.Eddy.VelErr = Quiver3VectorErr;
        handles.Eddy.VelSTDE = VelSTDE;
        handles.Eddy.VelSTDN = VelSTDN;
        handles.Eddy.VelSTDUp = VelSTDUp;
        handles.Eddy.VelSTDD = VelSTDD;
%         
        handles.Eddy.track = track(2:end)';
        handles.Eddy.depth = depth;
        handles.Eddy.Lat = Lat;
        handles.Eddy.Long = Long;   
        Result = EddyVisco(handles);
        LatEddy = [LatEddy Result.Lat];
        LongEddy = [LongEddy Result.Long];
        MUEddy = [MUEddy Result.MU];
        if handles.Eddy.Method ~= 2            
            VectorEddy = [VectorEddy Result.Vector];
            LatSimulEddy = [LatSimulEddy Result.LatitudeSimul];
            LongSimulEddy = [LongSimulEddy Result.LongitudeSimul];
            DepthSimulEddy = [DepthSimulEddy Result.DepthSimul];
        end
    end
        
%         EddyVisco(handles)
end
if handles.EddySimulV ==1 || handles.EddyExpV ==1
    VectorEddy.VectorEddy = VectorEddy;
    VectorEddy.LatSimulEddy=LatSimulEddy;
    VectorEddy.LongSimulEddy=LongSimulEddy;
    VectorEddy.DepthSimulEddy=DepthSimulEddy;
end
if handles.plotpreview == 1
    figure(plotprev);
    legend(pl,legendstr)
    xlabel('Longitude(deg)')
    ylabel('Latitude(deg)')
    title('Plot Preview')
end
% save('C:\Users\Roberta\Desktop\Ricardo\Variables Save\excelall.mat','excel')
QTotal = AreaTot.*TotalAverageVelocity'; 
% if PlSav == 1
%     close(gcf)
% end

if batimetriav == 0
%     assignin('base', 'LatTot', LatTot)
%     assignin('base', 'LongTot', LongTot)
%     assignin('base', 'depth1', depth1)
    NumofPointsBat = get(handles.Options3D.NumofPointsBat,'UserData');
    Xr = linspace(min(LatTot),max(LatTot),NumofPointsBat);
    Yr = linspace(min(LongTot),max(LongTot),NumofPointsBat);
    [XrR,YrR] = meshgrid(Xr,Yr);
    FrR = griddata(LongTot,LatTot,depth1,YrR,XrR);
else
    [Xr,Yr,Fr] = batimetriaToc(min(LatTot),min(LongTot),max(LatTot),max(LongTot),filenamebat, pathbat,handles);
    NumofPointsBat = get(handles.Options3D.NumofPointsBat,'UserData');
    [XrR,YrR] = meshgrid(linspace(min(LatTot),max(LatTot),NumofPointsBat),...
        linspace(min(LongTot),max(LongTot),NumofPointsBat));    
    FrR = griddata(Xr,Yr,Fr,XrR,YrR);
    FrR = FrR - A.GPS.Altitude(1);
end
if super == 1
    htest = figure;
    clf
    ax(3) = axes('Position',[0.25,0.1,0.7,0.7]);
    plot_google_map
    set(ax(3),'ytick',[]);set(ax(3),'xtick',[])
    ax(1) = axes('Position',[0.25,0.1,0.7,0.7]);
    [LowerLimQuiver, UpperLimQuiver,sQuiver] = limites(sqrt(AverageVelocityVectorTot1.^2 + AverageVelocityVectorTot2.^2));
    hold on
    xlabel('Longitude(deg)');ylabel('Latitude(deg)');

    if batimetriav == 0
%         NumofPointsBat = get(handles.Options3D.NumofPointsBat,'UserData');
%         Xr = linspace(min(LatTot),max(LatTot),NumofPointsBat);
%         Yr = linspace(min(LongTot),max(LongTot),NumofPointsBat);
%         [XrR,YrR] = meshgrid(Xr,Yr);
%         FrR = griddata(LongTot,LatTot,depth1,YrR,XrR);
        set(ax(1), 'XAxisLocation','bottom',...
                 'YAxisLocation','left',...
                 'Color','none');
        [~,h2] = contourf(YrR,XrR,FrR,get(handles.Options3D.ContourLines,'UserData'));       

        colormap(get(handles.Options3D.ContourMap,'UserData'))
        freezeColors        
        hold on
        for q=1:length(h2)
        set(h2(q),'LineStyle','none');
        end
        hleve=get(h2,'LevelList');
        hchil = get(h2,'Children');
        cont_level = zeros(1,length(hchil));
        legend_entries = cell(1,length(hleve)-1);
        TFtot = true(1,length(hchil));
        for dd= 1:length(hchil)
            cont_level(dd) = get(hchil(dd),'userdata');
            TF = cont_level(dd)==cont_level;
            if sum(TF(:)) == 1
                TFtot(dd) = true;
            else
                TFtot(dd) = false;
            end
            if hleve(1)==cont_level(dd)
                TFtot(dd) = false;
            end
        end
        for kg = 1:length(hleve)-1
            legend_entries{kg} = [num2str(hleve(kg),3) 'm<x<' num2str(hleve(kg+1),3) 'm'];
        end
        [~,I]=sort(cont_level(TFtot));
        B1 = hchil(TFtot);
        legend(B1(I),legend_entries,'Location','southwest')
        colormap(C)        
        if quiv3 == 1
            htest2 = figure;
            hold on
            scale= get(handles.Options3D.quiver3doption,'UserData');
            quiver3(Quiver3VectorXTot,Quiver3VectorYTot,Quiver3VectorZTot,Quiver3VectorUTot,Quiver3VectorVTot,Quiver3VectorWTot,scale)

            ylabel('Latitude (deg)')
            xlabel('Longitude (deg)')
            zlabel('Depth (m)')
            handle = surf(YrR,XrR,FrR);
            view(3)
            xlim([min(YrR(:)) max(YrR(:))])
            ylim([min(XrR(:)) max(XrR(:))])
            set(handle,'EdgeColor','none')
            colormap('gray')
        end
    else
%         [Xr,Yr,Fr] = batimetriaToc(min(LatTot),min(LongTot),max(LatTot),max(LongTot),filenamebat, pathbat,handles);
%         NumofPointsBat = get(handles.Options3D.NumofPointsBat,'UserData');
%         [XrR,YrR] = meshgrid(linspace(min(LatTot),max(LatTot),NumofPointsBat),...
%             linspace(min(LongTot),max(LongTot),NumofPointsBat));    
%         FrR = griddata(Xr,Yr,Fr,XrR,YrR);
%         FrR = FrR - A.GPS.Altitude(1);
        [~,hh]=contourf(YrR,XrR,FrR,get(handles.Options3D.ContourLines,'UserData'));
        for q=1:length(hh)
            set(hh(q),'LineStyle','none');
        end  
        colormap(get(handles.Options3D.ContourMap,'UserData'))
        freezeColors
        hleve=get(hh,'LevelList');
        hchil = get(hh,'Children');
        cont_level = zeros(1,length(hchil));
        legend_entries = cell(1,length(hleve)-1);
        TFtot = true(1,length(hchil));
        for dd= 1:length(hchil)
            cont_level(dd) = get(hchil(dd),'userdata');
            TF = cont_level(dd)==cont_level;
            if sum(TF(:)) == 1
                TFtot(dd) = true;
            else
                TFtot(dd) = false;
            end
            if hleve(1)==cont_level(dd)
                TFtot(dd) = false;
            end
        end
        for kg = 1:length(hleve)-1
            legend_entries{kg} = [num2str(hleve(kg),3) 'm<x<' num2str(hleve(kg+1),3) 'm'];
        end
        [~,I]=sort(cont_level(TFtot));
        B1 = hchil(TFtot);
        legend(B1(I),legend_entries,'Location','southwest')
        hold on
        set(ax(1), 'XAxisLocation','bottom',...
                     'YAxisLocation','left','Color','none');
        colormap(C)
        if quiv3 == 1
            htest2 = figure;
            scale= get(handles.Options3D.quiver3doption,'UserData');
            quiver3(Quiver3VectorXTot,Quiver3VectorYTot,Quiver3VectorZTot,Quiver3VectorUTot,Quiver3VectorVTot,Quiver3VectorWTot,scale,'Color','g')
            hold on
            ylabel('Latitude (deg)')
            xlabel('Longitude (deg)')
            zlabel('Depth (m)')
            handle = surf(YrR,XrR,FrR);
            view(3)
            xlim([min(YrR(:)) max(YrR(:))])
            ylim([min(XrR(:)) max(XrR(:))])
            set(handle,'EdgeColor','none')        
        end

    end
    if quiv3 == 1
        savefig(htest2,[Directory 'Analysis 3D\QuiverAll3D.fig']);
        clf
        close(htest2)
    end
    figure(htest)
    hold on    
    ax(2) = axes('Position',[0.25,0.1,0.7,0.7]);
    quiverwcolorbar(LongTot,LatTot,AverageVelocityVectorTot1,AverageVelocityVectorTot2,...
        get(handles.Options3D.VisTopQuiv,'UserData'),'bounds',[LowerLimQuiver, UpperLimQuiver]);
    set(gca,'Ydir','normal')

    
    h_bar = get(htest,'Children');
    h_bar=h_bar(1);
    set(h_bar,'Location','north')
    set(ax(2),'Color','none')
    set(ax(2),'ytick',[]);set(ax(2),'xtick',[])
%     h_bar = colorbar('North');
%     colormap(jet(64))
    set(get(h_bar,'xlabel'),'String', 'Velocity(m/s)');
    initpos = get(h_bar,'Position');
    initfontsize = get(h_bar,'FontSize');
    set(h_bar,'Position',[initpos(1)+initpos(3)/2 initpos(2) initpos(3)/2 initpos(4)/2],...
        'FontSize',initfontsize*.75)
%     caxis([LowerLimQuiver, UpperLimQuiver]);  

    linkaxes([ax(1),ax(2),ax(3)],'xy');
    ylim(ax(2),'auto')
    Ylim=get(ax(1),'Ylim');Xlim=get(ax(1),'Xlim');
%     save('C:\Users\Roberta\Desktop\Ricardo\Variables Save\var.mat','LatTot','LongTot','Ylim','Xlim','XrR','YrR')
    ax(4) = axes('Position',[0.1,0.75,0.23,0.23]);
    box on
    set(ax(4),'ytick',[]);set(ax(4),'xtick',[])
    limitsmini = get(handles.Options3D.smallgraphlimts,'UserData');
    ylim(ax(4),limitsmini(3:4));xlim(ax(4),limitsmini(1:2))
    hold on
    plot_google_map
    patch([Xlim(1) Xlim(1) Xlim(2) Xlim(2)],[Ylim(1) Ylim(2) Ylim(1) Ylim(2)],zeros(4,1),'FaceColor','none','Marker','o')
%     savefig(figure(htest),[Directory 'Analysis 3D\All Sections.fig']);
    
    if grid == 1 
%         h = figure;
        [Vr1,Vr2,Vr3,Lat1,Lon1,utmx, utmy,zone,TFTot] = LoadDataToGrid(LongTot,LatTot,zone,...
            AverageVelocityVectorTot1,AverageVelocityVectorTot2,...
            AverageVelocityVectorTot3,filenamegrid, pathgrid);
        idx = ~isnan(Vr1) & ~isnan(Vr2); 
        GridVar.Longitudeidx = Lon1(idx);
        GridVar.Latitudeidx = Lat1(idx);
        GridVar.utmxidx = utmx(idx);
        GridVar.utmyidx = utmy(idx);
        GridVar.VelocityEidx = Vr1(idx);
        GridVar.VelocityNidx = Vr2(idx);
        GridVar.VelocityUidx = Vr3(idx);
        GridVar.Longitude = Lon1;
        GridVar.Latitude = Lat1;
        GridVar.utmx = utmx;
        GridVar.utmy = utmy;
        GridVar.zone = zone;
        GridVar.VelocityE = Vr1;
        GridVar.VelocityN = Vr2;
        GridVar.VelocityU = Vr3;
        GridVar.idx = and(TFTot(:),idx(:));
        
        save([Directory 'Analysis 3D\GridData.mat'],'GridVar')
        
        close(htest)
        htest = figure;
        ax1(3) = axes('Position',[0.25,0.1,0.7,0.7]);
        plot_google_map
        set(ax1(3),'ytick',[]);set(ax1(3),'xtick',[])
        ax1(1) = axes('Position',[0.25,0.1,0.7,0.7]);
        set(ax1(1), 'XAxisLocation','bottom',...
                 'YAxisLocation','left','Color','none');

        hold on
        xlabel('Longitude(deg)');ylabel('Latitude(deg)');
%         if batimetriav ==1
        [~,hh]=contourf(YrR,XrR,FrR,'LineStyle','none');
        for q=1:length(hh)
            set(hh(q),'LineStyle','none');
        end  
        colormap(get(handles.Options3D.ContourMap,'UserData'))
        freezeColors
        hleve=get(hh,'LevelList');
        hchil = get(hh,'Children');
        cont_level = zeros(1,length(hchil));
        legend_entries = cell(1,length(hleve)-1);
        TFtot = true(1,length(hchil));
        for dd= 1:length(hchil)
            cont_level(dd) = get(hchil(dd),'userdata');
            TF = cont_level(dd)==cont_level;
            if sum(TF(:)) == 1
                TFtot(dd) = true;
            else
                TFtot(dd) = false;
            end
            if hleve(1)==cont_level(dd)
                TFtot(dd) = false;
            end
        end
        for kg = 1:length(hleve)-1
            legend_entries{kg} = [num2str(hleve(kg),3) 'm<x<' num2str(hleve(kg+1),3) 'm'];
        end
        [~,I]=sort(cont_level(TFtot));
        B1 = hchil(TFtot);
        legend(B1(I),legend_entries,'Location','southwest')
        hold on
        hold on

        ax1(2) = axes('Position',[0.25,0.1,0.7,0.7]);
        quiverwcolorbar(Lon1(idx),Lat1(idx),Vr1(idx),Vr2(idx),...
            get(handles.Options3D.GridTopQuiv,'UserData'),'bounds',[LowerLimQuiver, UpperLimQuiver]);
        set(gca,'Ydir','normal')
        
        h_bar1 = get(htest,'Children');
        h_bar1=h_bar1(1);
        set(h_bar1,'Location','north')
        set(ax1(2),'Color','none')
        set(ax1(2),'ytick',[]);set(ax1(2),'xtick',[])
        colormap(jet(64))
        set(get(h_bar1,'xlabel'),'String', 'Velocity(m/s)');
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
        limitsmini = get(handles.Options3D.smallgraphlimts,'UserData');
        ylim(ax1(4),limitsmini(3:4));xlim(ax1(4),limitsmini(1:2))
        hold on
        plot_google_map
        patch([Xlim(1) Xlim(1) Xlim(2) Xlim(2)],[Ylim(1) Ylim(2) Ylim(1) Ylim(2)],zeros(4,1),'FaceColor','none','Marker','o')
        savefig(figure(htest),[Directory 'Analysis 3D\All Section Grid.fig']);
    end


else
    if quiv3 ==1 && batimetriav ==1
%         numbat = get(handles.Options3D.NumofPointsBat,'UserData');
%         [Xr,Yr,Fr] = batimetriaToc(min(LatTot),min(LongTot),max(LatTot),max(LongTot),filenamebat, pathbat,handles);
%         [XrR,YrR] = meshgrid(linspace(min(LatTot),max(LatTot),numbat),...
%             linspace(min(LongTot),max(LongTot),numbat));    
%         FrR = griddata(Xr,Yr,Fr,XrR,YrR);
%         FrR = FrR - A.GPS.Altitude(1);
        htest2 = figure;
        scale= get(handles.Options3D.quiver3doption,'UserData');
        quiver3(Quiver3VectorXTot,Quiver3VectorYTot,Quiver3VectorZTot,Quiver3VectorUTot,Quiver3VectorVTot,Quiver3VectorWTot,scale,'Color','g')
        hold on
        ylabel('Latitude (deg)')
        xlabel('Longitude (deg)')
        zlabel('Depth (m)')
        handle = surf(YrR,XrR,FrR);
        set(handle,'EdgeColor','none')
        savefig(htest2,[Directory 'Analysis 3D\QuiverAll3D.fig']);
        
    elseif quiv3 ==1 && batimetriav ==0
%         numbat = get(handles.Options3D.NumofPointsBat,'UserData');
%         Xr = linspace(min(LatTot),max(LatTot),numbat);
%         Yr = linspace(min(LongTot),max(LongTot),numbat);
%         [XrR,YrR] = meshgrid(Xr,Yr);
%         FrR = griddata(LongTot,LatTot,depth1,YrR,XrR);
        htest2 = figure;
        hold on
        scale= get(handles.Options3D.quiver3doption,'UserData');
        quiver3(Quiver3VectorXTot,Quiver3VectorYTot,Quiver3VectorZTot,Quiver3VectorUTot,Quiver3VectorVTot,Quiver3VectorWTot,scale)

        ylabel('Latitude (deg)')
        xlabel('Longitude (deg)')
        zlabel('Depth (m)')
        handle = surf(YrR,XrR,FrR);
        set(handle,'EdgeColor','none')
        colormap;
        savefig(htest2,[Directory 'Analysis 3D\QuiverAll3D.fig']);
        
    end
    if simul == 1 && batimetriav ==1
%         numbat = get(handles.Options3D.NumofPointsBat,'UserData');
%         [Xr,Yr,Fr] = batimetriaToc(min(LatTot),min(LongTot),max(LatTot),max(LongTot),filenamebat, pathbat,handles);
%         [XrR,YrR] = meshgrid(linspace(min(LatTot),max(LatTot),numbat),...
%             linspace(min(LongTot),max(LongTot),numbat));    
%         FrR = griddata(Xr,Yr,Fr,XrR,YrR);
%         FrR = FrR - A.GPS.Altitude(1);
    elseif simul == 1 && batimetriav ==0
%         numbat = get(handles.Options3D.NumofPointsBat,'UserData');
%         Xr = linspace(min(LatTot),max(LatTot),numbat);
%         Yr = linspace(min(LongTot),max(LongTot),numbat);
%         [XrR,YrR] = meshgrid(Xr,Yr);
%         FrR = griddata(LongTot,LatTot,depth1,YrR,XrR);
        IN = inpolygon(Xr,Yr,[min(LongTot) max(LongTot)],[min(LatTot) max(LatTot)]);
        FrR(~IN) = NaN;
    elseif handles.ModelComparisonV == 1 && grid == 1
            numbat = get(handles.Options3D.NumofPointsBat,'UserData');
            Xr = linspace(min(LatTot),max(LatTot),numbat);
            Yr = linspace(min(LongTot),max(LongTot),numbat);
            [XrR,YrR] = meshgrid(Xr,Yr);
            FrR = griddata(LongTot,LatTot,depth1,YrR,XrR);
            IN = inpolygon(Xr,Yr,[min(LongTot) max(LongTot)],[min(LatTot) max(LatTot)]);
            FrR(~IN) = NaN;
%     else
%         XrR = [];YrR=[];FrR=[];
    end
    if grid ==1 && StreamLines == 1
        [Vr1,Vr2,Vr3,Lat1,Lon1] = LoadDataToGrid(LongTot,LatTot,zone,...
            AverageVelocityVectorTot1,AverageVelocityVectorTot2,...
            AverageVelocityVectorTot3,filenamegrid, pathgrid);
        idx = ~isnan(Vr1) & ~isnan(Vr2); 
        GridVar.Longitudeidx = Lon1(idx);
        GridVar.Latitudeidx = Lat1(idx);
        GridVar.VelocityEidx = Vr1(idx);
        GridVar.VelocityNidx = Vr2(idx);
        GridVar.VelocityUidx = Vr3(idx);
        GridVar.Longitude = Lon1;
        GridVar.Latitude = Lat1;
        GridVar.VelocityE = Vr1;
        GridVar.VelocityN = Vr2;
        GridVar.VelocityU = Vr3;
        save([Directory 'Analysis 3D\GridData.mat'],'GridVar')
    elseif StreamLines ==1
        numbat = get(handles.Options3D.NumofPointsBat,'UserData');
        Xr = linspace(min(LatTot),max(LatTot),numbat);
        Yr = linspace(min(LongTot),max(LongTot),numbat);
        [XrR,YrR] = meshgrid(Xr,Yr);
        Vr1 = griddata(LongTot,LatTot,AverageVelocityVectorTot1,YrR,XrR);
        Vr2 = griddata(LongTot,LatTot,AverageVelocityVectorTot2,YrR,XrR);
        Vr3 = griddata(LongTot,LatTot,AverageVelocityVectorTot3,YrR,XrR);
        idx = ~isnan(Vr1) & ~isnan(Vr2); 
        GridVar.Longitudeidx = YrR(idx);
        GridVar.Latitudeidx = XrR(idx);
        GridVar.VelocityEidx = Vr1(idx);
        GridVar.VelocityNidx = Vr2(idx);
        GridVar.VelocityUidx = Vr3(idx);
        GridVar.Longitude = YrR;
        GridVar.Latitude = XrR;
        GridVar.VelocityE = Vr1;
        GridVar.VelocityN = Vr2;
        GridVar.VelocityU = Vr3;
        
    end
end
if kml3 == 1 
    display('Writting KML...')
    copyfile([pwd '\googleearth\data\redcone.dae'],[Directory 'Analysis 3D\redcone.dae'])
    GoogleEarthQuiver([Directory 'Analysis 3D\Todas Seçoes'],LongTot,LatTot,AverageVelocityVectorTot1,...
                AverageVelocityVectorTot2,AverageVelocityVectorTot3,handles); 
    display('Done!')
end
if handles.ModelComparisonV == 1 && grid == 1
    if super == 1 || StreamLines == 1
    else        
        [Vr1,Vr2,Vr3,Lat1,Lon1] = LoadDataToGrid(LongTot,LatTot,zone,...
            AverageVelocityVectorTot1,AverageVelocityVectorTot2,...
            AverageVelocityVectorTot3,filenamegrid, pathgrid);
        idx = ~isnan(Vr1) & ~isnan(Vr2); 
        GridVar.Longitudeidx = Lon1(idx);
        GridVar.Latitudeidx = Lat1(idx);
        GridVar.VelocityEidx = Vr1(idx);
        GridVar.VelocityNidx = Vr2(idx);
        GridVar.VelocityUidx = Vr3(idx);
        GridVar.Longitude = Lon1;
        GridVar.Latitude = Lat1;
        GridVar.VelocityE = Vr1;
        GridVar.VelocityN = Vr2;
        GridVar.VelocityU = Vr3;
        save([Directory 'Analysis 3D\GridData.mat'],'GridVar')
    end
    QQComp = load(handles.ModelComparisionVa.FilenamePath);
    XComp = QQComp.data.XComp;
    YComp = QQComp.data.YComp;
    htest = figure;
    ax1(1) = axes('Position',[0.25,0.1,0.7,0.7]);
    plot_google_map
    ax1(2) = axes('Position',[0.25,0.1,0.7,0.7]);
    XComp = XComp(idx);
    YComp = YComp(idx);
    magD= sqrt(Vr1(idx).^2+Vr2(idx).^2);
    MagModel = sqrt(XComp.^2+YComp.^2);
    if handles.ModelComparisionVa.Type == 1 %% Error
        Error = abs(magD-MagModel)./magD;
        scatter(Lon1(idx),Lat1(idx),[],Error,'filled');
        [~, UpperLim,s] = limites(Error);
        LowerLim=0;
        UpperLim = UpperLim-1.5*s;
    else %%Difference
        Difference = magD-MagModel;
        scatter(Lon1(idx),Lat1(idx),[],Difference,'filled');
        [LowerLim, UpperLim,~] = limites(Difference);
    end
    set(gca,'Ydir','normal')

    h_bar1 = colorbar('north');
    caxis([LowerLim, UpperLim])
    set(ax1(2),'Color','none')
    set(ax1(2),'ytick',[]);set(ax1(2),'xtick',[])
    colormap(jet(64))
    if handles.ModelComparisionVa.Type == 1
        set(get(h_bar1,'xlabel'),'String', 'Error(%)');
    else
        set(get(h_bar1,'xlabel'),'String', 'Velocity Difference(m/s)');
    end
    initpos2 = get(h_bar1,'Position');
    initfontsize2 = get(h_bar1,'FontSize');
    set(h_bar1,'Position',[initpos2(1)+initpos2(3)/2 initpos2(2) initpos2(3)/2 initpos2(4)/2],...
        'FontSize',initfontsize2*.75)

    linkaxes([ax1(1),ax1(2)],'xy');
    ylim(ax1(2),'auto')
    xlim(ax1(2),'auto')
    Ylim=get(ax1(2),'Ylim');Xlim=get(ax1(2),'Xlim');
    ax1(3) = axes('Position',[0.1,0.75,0.23,0.23]);
    box on
    set(ax1(3),'ytick',[]);set(ax1(3),'xtick',[])
    limitsmini = get(handles.Options3D.smallgraphlimts,'UserData');
    ylim(ax1(3),limitsmini(3:4));xlim(ax1(3),limitsmini(1:2))
    hold on
    plot_google_map
    patch([Xlim(1) Xlim(1) Xlim(2) Xlim(2)],[Ylim(1) Ylim(2) Ylim(1) Ylim(2)],...
        zeros(4,1),'FaceColor','none','Marker','o')     
end
if handles.Eddy3DV == 1
    if grid ==1
        figure
        scatter(LongEddy,LatEddy,[],MUEddy,'filled');
        [MuGrid,Lat1,Lon1,utmx, utmy,zone,TFTot] = LoadDataToGridEddy(LongEddy,LatEddy,zone,...
            MUEddy,filenamegrid, pathgrid);        
        hEddy = figure;
        axEddy(1) = axes('Position',[0.25,0.1,0.7,0.7]);
        plot_google_map
        axEddy(2) = axes('Position',[0.25,0.1,0.7,0.7]);
        
        [LowerLimMu, UpperLimMu,~] = limites(MuGrid);
        h=pcolor(Lon1,Lat1,MuGrid);
        
        
        assignin('base','Lon',Lon1)
        assignin('base','Lat',Lat1)
        assignin('base','MuGrid',MuGrid)
        
        
        xlabel('Longitude(deg)');ylabel('Latitude(deg)')
        set(h,'EdgeColor','none')
        set(gca,'Ydir','normal')
        h_barEddy = colorbar('north');
        set(axEddy(2),'Color','none')
        set(axEddy(2),'ytick',[]);set(axEddy(2),'xtick',[])
        set(get(h_barEddy,'xlabel'),'String', 'Eddy Viscosity {\nu}_t (m^2.s^{-1})');
        initpos2 = get(h_barEddy,'Position');
        initfontsize2 = get(h_barEddy,'FontSize');
        set(h_barEddy,'Position',[initpos2(1)+initpos2(3)/2 initpos2(2) initpos2(3)/2 initpos2(4)/2],...
            'FontSize',initfontsize2*.75)
        caxis([0, UpperLimMu])

        linkaxes([axEddy(1),axEddy(2)],'xy');
        ylim(axEddy(2),'auto')
        xlim(axEddy(2),'auto')
        Ylim=get(axEddy(2),'Ylim');Xlim=get(axEddy(2),'Xlim');
        axEddy(3) = axes('Position',[0.1,0.75,0.23,0.23]);
        box on
        set(axEddy(3),'ytick',[]);set(axEddy(3),'xtick',[])
        limitsmini = get(handles.Options3D.smallgraphlimts,'UserData');
        ylim(axEddy(3),limitsmini(3:4));xlim(axEddy(3),limitsmini(1:2))
        hold on
        plot_google_map
        patch([Xlim(1) Xlim(1) Xlim(2) Xlim(2)],[Ylim(1) Ylim(2) Ylim(1) Ylim(2)],...
            zeros(4,1),'FaceColor','none','Marker','o')     
      


        
    else
        
        
    end
    
    
end
alt = A.GPS.Altitude(1);
toc





    




