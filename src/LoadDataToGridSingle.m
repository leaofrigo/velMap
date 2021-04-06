function [Vr1,Vr2,Vr3,Lat1,Lon1,utmx, utmy,zone,TFTot] = LoadDataToGridSingle(LongTot,LatTot,zone,...
    AverageVelocityVectorTot1,AverageVelocityVectorTot2,AverageVelocityVectorTot3,filename, path)
% [filename, path]= uigetfile({'*.mat'},'Selecione um arquivo com o grid desejado');
QQ = load([path filename]);
x = QQ.data.X;
y = QQ.data.Y;
SizeMat = size(x);
if ~isempty(zone)
    Zone = repmat(zone,SizeMat(1),1);
    Lat1 = zeros(SizeMat);
    Lon1 = zeros(SizeMat);
    for i = 1:SizeMat(2)
        [LatTEMP,Lon1TEMP] = utm2deg(x(:,i),y(:,i),Zone);
        Lat1(:,i) = LatTEMP;
        Lon1(:,i) = Lon1TEMP;
    end
    utmx=x;
    utmy=y;
else
    Lat1 = x;
    Lon1 = y;
    utmx = zeros(SizeMat);
    utmy = zeros(SizeMat);
    zone = utmy;
    for j = 1:SizeMat(2)
        [utmxtemp,utmytemp,zonetemp] = deg2utm(Lat1(:,j),Lon1(:,j));
        utmx(:,j) = utmxtemp;
        utmy(:,j) = utmytemp;
        zone(:,j) = zonetemp;
    end
end
TFLat = and(Lat1>=min(min(LatTot)),Lat1<=max(max(LatTot)));
TFLon = and(Lon1>=min(min(LongTot)),Lon1<=max(max(LongTot)));
TFTot = and(TFLat==1,TFLon==1);
Lat1(~TFTot) = nan;
Lon1(~TFTot) = nan;

Vr1 = griddata(LongTot,LatTot,AverageVelocityVectorTot1,Lon1,Lat1,'nearest');
Vr2 = griddata(LongTot,LatTot,AverageVelocityVectorTot2,Lon1,Lat1,'nearest');
Vr3 = griddata(LongTot,LatTot,AverageVelocityVectorTot3,Lon1,Lat1,'nearest');