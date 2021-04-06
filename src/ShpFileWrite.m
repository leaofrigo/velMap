function [] = ShpFileWrite(Long,Lat,FileName,handles)

if handles.run3d == 1
    type = get(handles.Options3D.ExpShp3d,'UserData');
    output = get(handles.Options3D.shpfileexp3d,'UserData');
else
    type = get(handles.Options2D.ExpShp2d,'UserData');
    output = get(handles.Options2D.shpfileexp,'UserData');
end
switch type
    case 'Geovector'
        P = geoshape(Lat,Long);
    case 'Geopoint'
        P = geopoint(Lat,Long);
end
shapewrite(P,FileName)
switch output
    case 'utm'
        formatSpec = '%6.1f, %7.0f; ';
        [x,y] = deg2utm(Lat,Long);
        xy = zeros(length(x)*2,1);
        xy(1:2:end-1) = x;
        xy(2:2:end)= y;
    case 'Lat & Long'
        formatSpec = '%3.7f, %3.7f; ';
        xy = zeros(length(Lat)*2,1);
        xy(1:2:end-1) = Lat;
        xy(2:2:end)= Long;
end



fid = fopen([FileName(1:end-4) '.txt'],'wt');
fprintf(fid,formatSpec,xy);
fclose(fid);

