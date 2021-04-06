function handles = PlotMultiSessionIn1(handles)
OpenFiles = handles.OpenFiles;
Directory = handles.Directory;
path = [Directory 'Multiplas Seções em uma Sessions'];
for rr =1:length(OpenFiles)
    path = [path '_' OpenFiles{rr}];
end
if length(path)>200
    path = path(1:200);
end
if exist(path,'dir')
else
    mkdir(path)
end
quiv2D = handles.Quiverplot.quiv2D;
replot = replot4gui(handles);
Long = replot.Long;
Lat = replot.Lat;
QuiverX = replot.Quiver3X;
QuiverY = replot.Quiver3Y;
QuiverU = replot.Quiver3U;
QuiverV = replot.Quiver3V;
StartEdge = replot.StartEdge;
Direction = handles.VelPopUp.V;
plotting = handles.plotsess.V;

Quiver3allX = replot.Quiver3allX;
Quiver3allY = replot.Quiver3allY;
Quiver3allZ = replot.Quiver3allZ;
Quiver3allU = replot.Quiver3allU;
Quiver3allV = replot.Quiver3allV;
Quiver3allW = replot.Quiver3allW;
Sumnum = replot.Sumnum;
r = zeros(length(Long),1);
mineachX = r; mineachY = r; maxeachX = r; maxeachY = r;
for q = 1:length(Long)
    if StartEdge{q} == 0
        mineachX(q) = QuiverX{q}(1);
        mineachY(q) = QuiverY{q}(1);
        maxeachX(q) = QuiverX{q}(end);
        maxeachY(q) = QuiverY{q}(end);
    else
        mineachX(q) = QuiverX{q}(end);
        mineachY(q) = QuiverY{q}(end);
        maxeachX(q) = QuiverX{q}(1);
        maxeachY(q) = QuiverY{q}(1);
    end
end
NumOfPoints = get(handles.options.NumPointsX,'UserData');
aveXmin = mean(mineachX);
aveYmin = mean(mineachY);
aveXmax = mean(maxeachX);
aveYmax = mean(maxeachY);
X1 = linspace(aveXmin,aveXmax,NumOfPoints);
Y1 = linspace(aveYmin,aveYmax,NumOfPoints);
TotalDist = deg2distance([aveYmin aveYmax],[aveXmin aveXmax]);



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
SumNumTot =[];output= cell(length(OpenFiles),1);outputQuiv1 = cell(length(OpenFiles),1);
outputQuiv2 = cell(length(OpenFiles),1);
outputKMLN = cell(length(OpenFiles),1); outputKMLE = cell(length(OpenFiles),1);
outputKMLU = cell(length(OpenFiles),1);DepthAll=cell(length(OpenFiles),1);
TrackAll=cell(length(OpenFiles),1);

NumSmooth = get(handles.options.NumPointsX,'UserData');
NumOfCellsY = NumSmooth/2;
NumOfCellsX = NumSmooth;
for n=1:length(OpenFiles)
    A = load ([Directory OpenFiles{n}]);
    velocidades = A.WaterTrack.Velocity;
    depth = - A.Summary.Depth;
    %Calcular distancia viajada pelo barco
    Dist = A.Summary.Track;
    y = diff([0 0; Dist]); % calcular diferenca entre pontos
    y1= sum(abs(y)); % somar a distancia absoluta entre pontos
    largura=norm(y1); % Achar a magnitude
    clear('velE','velN','velU','velD');
    velE(:,:) = velocidades(:,1,:);
    velN(:,:) = velocidades(:,2,:);
    velU(:,:) = velocidades(:,3,:);
    velD(:,:) = velocidades(:,4,:);
    velMag = sqrt(velE.^2+velN.^2+velU.^2);
    TamanhoMatrix = size(velocidades);
    StartEdge1 = A.Setup.startEdge;
    xx = zeros(1,TamanhoMatrix(3)+1);
    NumpFile(n+1) = NumpFile(n)+TamanhoMatrix(3);

    for i = 1:TamanhoMatrix(3)
        xx(i+1) = xx(i) + norm(y(i,:)); 
    end
    TFvel = ~isnan(velD);
    %Calculando um linha reta entre primeiro e ultimo ponto
    DistMag = xx(end) - xx(1);
    xxNorm = xx/DistMag;      
%     DistMag = sqrt((Dist(1,1)-Dist(end,1))^2+(Dist(1,2)-Dist(end,2))^2);
    xx = xxNorm*TotalDist;
    xxCentro = zeros(1,length(xx)-1);  
    ComecoCell = A.System.Cell_Start;
    CellSize = A.System.Cell_Size;
    NumOfCells  = A.Summary.Cells;


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
    %Go Through each column of data
    VelocityTopExtrp = zeros(ExtrapStruct.NumCellsExtrap,TamanhoMatrix(3));
    VelocityBotExtrp =VelocityTopExtrp; Ytop =VelocityTopExtrp; Ybot=VelocityTopExtrp;
    for k = 1:TamanhoMatrix(3)
        [Sumnum(k+1),Sumnum2(k+1),X,Y,AverageVelocityVector(k,:),AverageVelocityMag(k),...
            cellDepths,zz(k+1),cellDepthsCentro,...
            ~,~,~,~,VelMagTemp,...
            Quiver3VectorX,Quiver3VectorY,Quiver3VectorZ,~] =...
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
%     [Face, VerticeXXYY,XXQuiver,YYQuiver] = CriarImagem(cellDepths,track,Sumnum,Sumnum2);
    [Area0,depthnew,tracknew] = edges(EdgeDist0,depth,type0,track,0,handles);
    [Area1,depthnew1,tracknew1] = edges(EdgeDist1,depthnew,type1,tracknew,1,handles);  
    DepthAll{n} = depth;
    TrackAll{n} = track;
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
%     
%     
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
    switch Direction
        case 'North'
            VelVec = Quiver3VectorVel(:,1);
        case 'East'
            VelVec = Quiver3VectorVel(:,2);
        case 'Upstream'
            VelVec = Quiver3VectorVel(:,3);
        case 'Transverse'
            VelVec = VelRL(:,1);
        case 'Longitudinal'
            VelVec = VelRL(:,2);
        case 'Magnitude'
            VelVec = Quiver3VectorVel(:,5);
        case 'Error'
            VelVec = Quiver3VectorVel(:,4);
    end

    output{n} = data2grid4multisess(xxCentro,NumOfCellsX,cellDepths,cellDepthsCentro,NumOfCellsY,...
        VelVec,Sumnum2,xx,tracknew1,depthnew1);
    if quiv2D{1} ==1
        
        switch quiv2D{2}
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
        outputQuiv1{n} = data2grid4multisess(xxCentro,NumOfCellsX,cellDepths,cellDepthsCentro,NumOfCellsY,...
            Vel1,Sumnum2,xx,tracknew1,depthnew1);
        switch quiv2D{3}
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
        outputQuiv2{n} = data2grid4multisess(xxCentro,NumOfCellsX,cellDepths,cellDepthsCentro,NumOfCellsY,...
            Vel2,Sumnum2,xx,tracknew1,depthnew1);          
    end
    if handles.kmlexp.V == 1 || handles.GridExp.V==1           
        outputKMLN{n} = data2grid4multisess(xxCentro,NumOfCellsX,cellDepths,cellDepthsCentro,...
            NumOfCellsY,Quiver3VectorV,Sumnum2,xx,tracknew1,depthnew1); 
        outputKMLE{n} = data2grid4multisess(xxCentro,NumOfCellsX,cellDepths,cellDepthsCentro,...
            NumOfCellsY,Quiver3VectorU,Sumnum2,xx,tracknew1,depthnew1);
        outputKMLU{n} = data2grid4multisess(xxCentro,NumOfCellsX,cellDepths,cellDepthsCentro,...
            NumOfCellsY,Quiver3VectorW,Sumnum2,xx,tracknew1,depthnew1);
    end
    AreaBody = abs(FindCrossAreaRiver(depth,xx));
    AreaTot(n) = Area0+Area1+AreaBody;
    Perimeter(n) = FindPerimeterRiver(depthnew1,tracknew1);
    HidrRad(n) = AreaTot(n)/Perimeter(n);  
    if handles.excelgen.V == 1
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
    end

end  
% save('C:\Users\Roberta\Desktop\Ricardo\Variables Save\Depths','DepthAll','TrackAll')
xxlinsp = linspace(min(xx),max(xx),NumOfCellsX);

% xr = linspace(min(x),max(x),round(length(x)/2));
MinDepth=1000;
for k = 1:length(output)
    wq = output{k,1}.depthmin;
    if wq<MinDepth
        MinDepth = wq;
    end
end
yylinsp = linspace(0,MinDepth,NumOfCellsY);
[xxgrid,yygrid] = meshgrid(xxlinsp,yylinsp);
Size = size(xxgrid);
VELTOT = zeros(Size(1),Size(2),length(output));
if quiv2D{1} ==1
    VELTOTQuiv1 = VELTOT;
    VELTOTQuiv2 = VELTOT;
end
Track = [];
Depth = [];
for k = 1:length(output)
    xtemp = output{k,1}.xxCentro;
    ytemp = output{k,1}.cellDepthsGrid;
    VelTemp = output{k,1}.AAgridVel;
    TF = and(~isnan(xtemp),~isnan(ytemp));
    if StartEdge{k} == 1
        xtemp = abs(xtemp-max(xtemp(:)));
        TrackAll{k} = abs(TrackAll{k}-max(TrackAll{k}));
    end
    Track = [Track TrackAll{k}(2:end)];
    Depth = [Depth DepthAll{k}'];
    if quiv2D{1} ==1
        Vel1Temp = outputQuiv1{k,1}.AAgridVel;
        Vel2Temp = outputQuiv2{k,1}.AAgridVel;
        VELTOTQuiv1(:,:,k) = griddata(xtemp(TF),ytemp(TF),...
            Vel1Temp(TF),xxgrid,yygrid);
        VELTOTQuiv2(:,:,k) = griddata(xtemp(TF),ytemp(TF),...
            Vel2Temp(TF),xxgrid,yygrid);        
    end
    if handles.kmlexp.V == 1 || handles.GridExp.V==1
        KMLTempN = outputKMLN{k,1}.AAgridVel;
        KMLTempE = outputKMLE{k,1}.AAgridVel;
        KMLTempU = outputKMLU{k,1}.AAgridVel;
        VELTOTKLMN(:,:,k) = griddata(xtemp(TF),ytemp(TF),...
            KMLTempN(TF),xxgrid,yygrid);
        VELTOTKLME(:,:,k) = griddata(xtemp(TF),ytemp(TF),...
            KMLTempE(TF),xxgrid,yygrid); 
        VELTOTKLMU(:,:,k) = griddata(xtemp(TF),ytemp(TF),...
            KMLTempU(TF),xxgrid,yygrid); 
    end
    VELTOT(:,:,k) = griddata(xtemp(TF),ytemp(TF),...
        VelTemp(TF),xxgrid,yygrid);
end
xr = linspace(min(Track),max(Track),round(length(Track)/k));
% save('C:\Users\Roberta\Desktop\Ricardo\Variables Save\track','Track','Depth')
% [Track,ia,~] = unique(Track);
% Depth =Depth(ia);
[Track,Depth] = consolidator(Track,Depth,@nanmean);
pp= spline(Track,Depth,xr);
pps = smooth(pp);
ppss = smooth(pps);
[~,depthnew,tracknew] = edges1(EdgeDist0,ppss,type0,xr,0,handles);
[~,depthnew1,tracknew1] = edges1(EdgeDist1,depthnew,type1,tracknew,1,handles);
VELTOT1 = nanmean(VELTOT,3);
if plotting == 1
    figu = figure;
    pp = pcolor(xxgrid,yygrid,VELTOT1);
    set(pp,'EdgeColor','none')
    [LowerVN,UpperVN] = limites(VELTOT1);
    colorbar('south')
    C = hsv(250);
    colormap(flipud(C([1:4:16,14:50,50:2:end-125,end-125:end-50,end-50:4:end-30],:))) % make the color bar look like the River Servoyor Live
    h_bar = findobj(gcf,'Tag','Colorbar');
    set(get(h_bar,'xlabel'),'String', 'Velocity(m/s)');
    initpos = get(h_bar,'Position');
    initfontsize = get(h_bar,'FontSize');
    set(h_bar,'Position',[initpos(1)+initpos(3)/2 initpos(2) initpos(3)/2 initpos(4)/2],...
        'FontSize',initfontsize*.75)
    hold on
    caxis([LowerVN UpperVN]);
    xlabel('Length (m)')
    ylabel('Depth (m)')
    hold on
    plot(tracknew1,depthnew1);

    set(figu,'Position',get(handles.Options2D.plotpostion,'UserData'))
    xlim([min(tracknew1) max(tracknew1)])
    ylim([min(depthnew1)+min(depthnew1)*.1 0])
   
    switch Direction
        case 'North'
            title('Velocity North');
            savefig(gcf,[path '\Velocity_North.fig']);
        case 'East'
            title('Velocity East');
            savefig(gcf,[path '\Velocity_East.fig']);
        case 'Upstream'
            title('Velocity Upstream');
            savefig(gcf,[path '\Velocity_Upstream.fig']);
        case 'Error'
            title('Error Velocity');
            savefig(gcf,[path '\Velocity_Error.fig']);
        case 'Magnitude'
            title('Magnitude da velocidade');
            savefig(gcf,[path '\Velocity_Magnitude.fig']);
        case 'Transverse'
            title('Velocity Transverse')
            savefig(gcf,[path '\Velocity_Transverse.fig']);
        case 'Longitudinal'
            title('Velocity Longitudinal')
            savefig(gcf,[path '\Velocity_Longitudinal.fig']);
    end
    if quiv2D{1} ==1
        VELTOTQuiv1 = nanmean(VELTOTQuiv1,3);
        VELTOTQuiv2 = nanmean(VELTOTQuiv2,3);
        hold on
        QuivMulti = get(handles.Options3D.QuivMulti,'UserData');
        quiver(xxgrid,yygrid,VELTOTQuiv1,VELTOTQuiv2,QuivMulti)
    end
end
if handles.excelgen.V == 1
    excel(n+1).section = {'Average All'};
    excel(n+1).Total_Q = mean([excel(1:n).Total_Q]);    
    excel(n+1).Ave_speed = mean([excel(1:n).Ave_speed]); 
    excel(n+1).totaldist = mean([excel(1:n).totaldist]);
    excel(n+1).diffxx1xxend = mean([excel(1:n).diffxx1xxend]);
    excel(n+1).maxdepthstart = mean([excel(1:n).maxdepthstart]);
    excel(n+1).perimeter = mean([excel(1:n).perimeter]);
    excel(n+1).AreaTot = mean([excel(1:n).AreaTot]);
    excel(n+1).hidrrad = mean([excel(1:n).hidrrad]);
    excel(n+1).AreaMeas = mean([excel(1:n).AreaMeas]);
    [path '\ExcelReport.xls']
    writeExcelReport(excel,[path '\ExcelReport.xls'])
    
end
if handles.kmlexp.V == 1 || handles.GridExp.V==1
    VELTOTKLMN = nanmean(VELTOTKLMN,3);
    VELTOTKLME = nanmean(VELTOTKLME,3);
    VELTOTKLMU = nanmean(VELTOTKLMU,3);
    VELTOTKLMN = nanmean(VELTOTKLMN);
    VELTOTKLME = nanmean(VELTOTKLME);
    VELTOTKLMU = nanmean(VELTOTKLMU);
end
if handles.kmlexp.V == 1 
    display('Writting KML...')
    copyfile('D:\Users\Desktop\modelo ricardo VMT\Matlab scripts\googleearth\data\redcone.dae',[path '\redcone.dae'])
    GoogleEarthQuiver([path '\Multi Session'],X1,Y1,VELTOTKLME,...
                VELTOTKLMN,VELTOTKLMU,handles); 
    display('Done!')
end
if handles.ShpExp.V==1
    display('Writting Shape File...')
    FileName = [path '\Multi Seção'];
    ShpFileWrite(X1,Y1,FileName,handles)
    display('Shape File Done!')
end

if handles.GridExp.V==1
    display('Writting Grid File...')
    zone = handles.GridExp.zone;
    pathgrid = handles.GridExp.path;
    filenamegrid = handles.GridExp.filename;

    [Vr1,Vr2,Vr3,Lat1,Lon1,utmx, utmy,zone,TFTot] = LoadDataToGridSingle(X1,Y1,zone,...
        VELTOTKLME,VELTOTKLMN,VELTOTKLMU,filenamegrid, pathgrid);    
    idx = ~isnan(Vr1) & ~isnan(Vr2);
    Longitudeidx = Lon1(idx);
    Latitudeidx = Lat1(idx);
    utmxidx = utmx(idx);
    utmyidx = utmy(idx);
    VelocityEidx = Vr1(idx);
    VelocityNidx = Vr2(idx);
    VelocityUidx = Vr3(idx);
    
    minX = ones(1,length(X1));minXI = minX;
    for i=1:length(X1)
        minXtemp = abs(X1(i) - Longitudeidx);
        minYtemp = abs(Y1(i) - Latitudeidx);
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
    GridVar2D.PrimeiroIdx = and(TFTot(:),idx(:));
    GridVar2D.SegundoIdx = minXI;
    
    save([path '\GridDataSec.mat'],'GridVar2D')
    display('Grid Done!')
end
if handles.VortLate.V ==1
    
    handles.SideVort.yygrid = yygrid; 
    handles.SideVort.xxgrid = xxgrid;     
    handles.SideVort.Vel = VELTOT1;
    handles.SideVort.track = tracknew1;
    handles.SideVort.depth = depthnew1;

end

