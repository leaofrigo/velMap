function Result = EddyVisco(handles)

data = get(handles.Options2D.VortSide,'UserData');
NUMBOX = data(1)+1;
perc = data(2);
Face = handles.Eddy.Face;
VerticeXXYY=handles.Eddy.VerticeXXYY;
% Vel = magnitude
Vel = handles.Eddy.Vel;
%VelE = East Velocity // VelN = North Velocity // VelU = Upstream Velocity
%// VelErr = Error Velocity //VelSTD = All directions STD dev. Vel
VelE = handles.Eddy.VelE;
VelN = handles.Eddy.VelN;
VelUp = handles.Eddy.VelUp;
VelErr = handles.Eddy.VelErr;
VelSTDE = handles.Eddy.VelSTDE;
VelSTDN = handles.Eddy.VelSTDN;
VelSTDUp = handles.Eddy.VelSTDUp;
VelSTDD = handles.Eddy.VelSTDD;
%
track = handles.Eddy.track;
depth = handles.Eddy.depth;
Lat = handles.Eddy.Lat;
Long = handles.Eddy.Long;
% StartEdge = handles.SideVort.StartEdge;
ft=fittype('a*(exp(x*b))');
% global qh qq qpatch
disttot = max(track)-min(track);
distperc = disttot*perc;
x = Face(:,1:2);%h.x=x;
xx = VerticeXXYY(x(:,1),1);%h.xx=xx;
xx2 = VerticeXXYY(x(:,2),1);%h.xx2=xx2;

y=Face(:,2:3);%h.y=y;
yy = VerticeXXYY(y(:,1),2);%h.yy=yy;
yy2 = VerticeXXYY(y(:,2),2);%h.yy2=yy2;
yyy = [yy yy2];%h.yyy=yyy;
yyymean=nanmean(yyy')';%h.yyymean=yyymean;
Face1 = Face(:,1);
Binnumber=1;
Kappa = .41;
ro = 1000;
C_mu = 0.09;
% TABLE = zeros(10,round((max(track)-distperc-min(track))/distperc)+1);
Latitude = zeros(1,round((max(track)-distperc-min(track))/distperc)+1);
Longitude = Latitude;
MU = Latitude;
Vector=[];
LatitudeSimul = [];
LongitudeSimul =[];
DepthSimul=[];
C0 = .2; % *Add capability in menu* usually between 0.19 and 0.21
for val = min(track):distperc:max(track)-distperc
    TFx = and(xx>val,xx2<val+distperc);%h.TFx=TFx;        
    yTF1 = yy(TFx);%h.yTF1=yTF1;
    yTF2 = yy2(TFx);%h.yTF2=yTF2;
    yyTF = [yTF1 yTF2];%h.yyTF=yyTF;
    xTF1 = xx(TFx);%h.xTF1=xTF1;
    xTF2 = xx2(TFx);%h.xTF2=xTF2;
    xxTF = [xTF1 xTF2];%h.xxTF=xxTF;
%     xyTF = [xxTF yyTF];%h.xyTF=xyTF;
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
    TfVelSTDE = Area;
    TfVelSTDN=Area;
    TfVelSTDUp=Area;
    Areasum = 0;
    Qcell = Area;
    QcellK = Area;
    QAveK= zeros(1,length(boxes)-1);
    Qave= zeros(1,length(boxes)-1);
    yloc =Qave;  
    Du=yloc;
    uprime1=yloc;
    uPrimeK = yloc;
    tau = yloc;
    for j = 1:length(boxes)-1
        if j ==length(boxes)-1
            TFbox(:,j) = and(yyymean(:,1)<=boxes(j),yyymean(:,1)>=boxes(j+1)); 
        else
            TFbox(:,j) = and(yyymean(:,1)<=boxes(j),yyymean(:,1)>boxes(j+1)); 
        end
        TFbox(:,j) = and(TFbox(:,j),TFx);
        TFface{j} = Face1(logical(TFbox(:,j)));
        TFVel{j} = Vel(logical(TFbox(:,j)));
%         Calculate "ktemp" from std Velocities
        TfVelSTDE{j} = VelSTDE(logical(TFbox(:,j)));
        TfVelSTDN{j} = VelSTDN(logical(TFbox(:,j)));
        TfVelSTDUp{j} = VelSTDUp(logical(TFbox(:,j)));
        ktemp = (TfVelSTDE{j}.^2 + TfVelSTDN{j}.^2 + TfVelSTDUp{j}.^2)/2;
%         
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
%         finish Calculating ktemp
        QcellK{j} = ktemp.*AreaTemp;
        QAveK(j) = nanmean(QcellK{j});
        uPrimeK(j) = QAveK(j)/boxarea; %this is k
        tau(j) = C0*ro.*uPrimeK(j);
%        
        
        
    end
    ubar=nanmean(uprime1);
    Du = zeros(1,length(uprime1)-1);
    Dy=Du;
    un = cell(1,length(uprime1)-1);
    for kk = 1:length(Du)
        Du(kk) = uprime1(kk)-uprime1(kk+1);
        Dy(kk) = yloc(kk)-yloc(kk+1);
        un{kk} = ['U_' num2str(kk)];
    end
    Mom = Du.*Dy;

    Depthmean = nanmean(depth(and(track>val,track<val+distperc)));
    DepthminMinusMax = abs(min(depth(and(track>val,track<val+distperc)))-max(depth(and(track>val,track<val+distperc))));
    DepthSTD = std(depth(and(track>val,track<val+distperc)));
    Latitude(Binnumber) = nanmean(Lat(and(track>val,track<val+distperc)));
    Longitude(Binnumber) = nanmean(Long(and(track>val,track<val+distperc)));
    trye = [yloc;uprime1];
    trye(:,length(trye)+1)=[min(depth(and(track>val,track<val+distperc)));0];
    trye(1,:)=trye(1,:)-min(trye(1,:));
    StructName = ['bin' num2str(Binnumber)];
    if handles.Eddy.Method ==1
        f = fit(trye(2,:)',trye(1,:)',ft,'StartPoint', [3, 30]);
        coeff = coeffvalues(f);
        u_tau = Kappa/coeff(2);
        mu_turb = u_tau*Kappa*linspace(0,max(trye(1,:)),10);
    
    
%     TABLE(1,Binnumber) = val;
%     TABLE(2,Binnumber) = val+distperc;
%     TABLE(3,Binnumber) = Depthmean;
%     TABLE(4,Binnumber) = DepthSTD;
%     TABLE(5,Binnumber) = DepthminMinusMax;
%     TABLE(6,Binnumber) = u_tau;
%     TABLE(7,Binnumber) = mean(mu_turb);
        
    
        Method.(StructName).start = val;
        Method.(StructName).end = val+distperc;
        Method.(StructName).Depthmean = Depthmean;
        Method.(StructName).DepthSTD = DepthSTD;
        Method.(StructName).DepthminMinusMax = DepthminMinusMax;
        Method.(StructName).u_tau = u_tau;
        Method.(StructName).mu_turb =  mu_turb;
        Method.(StructName).mu_turbMEAN =  mean(mu_turb);
        Vector = [Vector mu_turb];
        LatitudeSimul = [LatitudeSimul ones(1,length(mu_turb))*Latitude(Binnumber)];
        LongitudeSimul = [LongitudeSimul ones(1,length(mu_turb))*Longitude(Binnumber)];
        DepthSimul = [DepthSimul linspace(min(depth(and(track>val,track<val+distperc))),...
            max(depth(and(track>val,track<val+distperc))),length(mu_turb))];
    elseif handles.Eddy.Method ==2
    
        MAXDuDy = FindDuDy(uprime1,yloc);
        Method.(StructName).start = val;
        Method.(StructName).end = val+distperc;
        Method.(StructName).Depthmean = Depthmean;
        Method.(StructName).DepthSTD = DepthSTD;
        Method.(StructName).DepthminMinusMax = DepthminMinusMax;
        Method.(StructName).u_tau = ubar; %% Until GUSTAVO comes up with something
        Method.(StructName).MAXDuDy = MAXDuDy;
        Method.(StructName).mu_turbMEAN =  ubar^2/MAXDuDy;
    
%     TABLE(8,Binnumber) = ro*ubar^2/MAXDuDy;
    elseif handles.Eddy.Method ==3
    
        Method.(StructName).start = val;
        Method.(StructName).end = val+distperc;
        Method.(StructName).Depthmean = Depthmean;
        Method.(StructName).DepthSTD = DepthSTD;
        Method.(StructName).DepthminMinusMax = DepthminMinusMax;
        Method.(StructName).u_tau = ubar; %% Until GUSTAVO comes up with something
        Method.(StructName).mu_turbMEAN = mean(ubar^2./abs(Du./Dy));
        Method.(StructName).mu_turb =  ubar^2./abs(Du./Dy);
        Vector = [Vector ubar^2./abs(Du./Dy)];
        LatitudeSimul = [LatitudeSimul ones(1,length(Du))*Latitude(Binnumber)];
        LongitudeSimul = [LongitudeSimul ones(1,length(Du))*Longitude(Binnumber)];
        DepthSimul = [DepthSimul linspace(min(depth(and(track>val,track<val+distperc))),...
            max(depth(and(track>val,track<val+distperc))),length(Du))];
    elseif handles.Eddy.Method ==4
    
%     TABLE(9,Binnumber) = mean(ro*ubar^2./abs(Du./Dy));
        K  = 1.5*mean((Du./Dy).^2);
        u_tauM4 = C_mu^(1/4)*K^(.5);
        mu_turbM4 = C_mu^(.25)*K^(1.5)*Kappa*linspace(0,max(trye(1,:)),10);
        Method.(StructName).start = val;
        Method.(StructName).end = val+distperc;
        Method.(StructName).Depthmean = Depthmean;
        Method.(StructName).DepthSTD = DepthSTD;
        Method.(StructName).DepthminMinusMax = DepthminMinusMax;
        Method.(StructName).u_tau = u_tauM4;
        Method.(StructName).mu_turb =  mu_turbM4;
        Method.(StructName).mu_turbMEAN =  mean(mu_turbM4);
    
    else
        uPrimeKtest = zeros(1,length(uPrimeK)-1);
        for tt = 1 : length(uPrimeKtest)
            uPrimeKtest(tt) = (uPrimeK(tt)+uPrimeK(tt+1))/2;            
        end
        LDu=Du;
        LDy=Dy;
        LuP=size(uPrimeKtest);
        lC0=C0;
        vt = C0.*uPrimeKtest./abs(Du./Dy);
        
        Method.(StructName).start = val;
        Method.(StructName).end = val+distperc;
        Method.(StructName).Depthmean = Depthmean;
        Method.(StructName).DepthSTD = DepthSTD;
        Method.(StructName).DepthminMinusMax = DepthminMinusMax;
        Method.(StructName).u_tau = (Du./Dy);
        Method.(StructName).mu_turb =  vt;
        Method.(StructName).mu_turbMEAN =  mean(vt);
        Vector = [Vector C0.*uPrimeKtest./abs(Du./Dy)];
        LatitudeSimul = [LatitudeSimul ones(1,length(Du))*Latitude(Binnumber)];
        LongitudeSimul = [LongitudeSimul ones(1,length(Du))*Longitude(Binnumber)];
        DepthSimul = [DepthSimul linspace(min(depth(and(track>val,track<val+distperc))),...
            max(depth(and(track>val,track<val+distperc))),length(Du))];
    end
    MU(Binnumber) = Method.(StructName).mu_turbMEAN;
    
%     TABLE(10,Binnumber) = mean(mu_turbM4);
    
    Binnumber = Binnumber+1;
        
end
Result.Lat = Latitude;
Result.Long = Longitude;
Result.MU = MU;
if handles.Eddy.Method ~=2 && handles.run2d ~= 1;
    Result.Vector = Vector;
    Result.LatitudeSimul = LatitudeSimul;
    Result.LongitudeSimul = LongitudeSimul;
    Result.DepthSimul = DepthSimul;
end
% assignin('base', 'Method', Result)
% assignin('base', 'TABLE', TABLE)
% assignin('base', 'uprime1', uprime1)
% assignin('base', 'yloc', yloc)
% assignin('base', 'Method2', Method2)
% assignin('base', 'Method3', Method3)
% assignin('base', 'Method4', Method4)

