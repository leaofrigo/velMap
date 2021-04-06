function [Lat, Long, Alt] = batimetriaToc(MinLat,MinLong,MaxLat,MaxLong,filename, path,handles)
% [filename, path]= uigetfile({'*.xyz'},'Selecione um arquivo com batimetria do rio desejado');
type = handles.zonebathV;
B = load([path filename]);
if type ==2

    [a,b,zone] = deg2utm(MinLat,MinLong);
    [a1,b1] = deg2utm(MaxLat,MaxLong);
    TFB = and(B(:,2) < b1,B(:,2) > b);
    TFA = and(B(:,1) < a1,B(:,1) > a);
    TF = and(TFB==1,TFA==1);
    [Lat,Long]= utm2deg(B(TF,1),B(TF,2),repmat(zone,size(B(TF,2))));
    Alt = B(TF,3);
else
    [a,b,zone] = deg2utm(MinLat,MinLong);
    [a1,b1] = deg2utm(MaxLat,MaxLong);
    [A2,A1]=deg2utm(B(:,1),B(:,2));
    B = [A1 A2];
    TFB = and(B(:,2) < b1,B(:,2) > b);
    TFA = and(B(:,1) < a1,B(:,1) > a);
    TF = and(TFB==1,TFA==1);
    [Lat,Long]= utm2deg(B(TF,1),B(TF,2),repmat(zone,size(B(TF,2))));
    Alt = B(TF,3);
end