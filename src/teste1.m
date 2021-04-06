kmlFileName = 'Tentativa1.kml';
Long = A.GPS.Longitude(1:end-1); 
Lat = A.GPS.Latitude(1:end-1);
Alt = A.GPS.Altitude(1:end-1);
u = AverageVelocityVector(:,1);
v = AverageVelocityVector(:,2);
w = AverageVelocityVector(:,3);
N = length(Long);
s_red = '';
d = date;
dnum = datevec( d );
arrowStr = 'redcone.dae';

for n = 1:(N-1)
    
    dnum2 = dnum;
    dnum3 = dnum;
    dnum2(5) = dnum(5) + n;
    dnum3(5) = dnum(5) + n +1;
    dstr = datestr( dnum2, 'yyyy-mm-ddTHH:MM:SSZ');
    dstr2 = datestr( dnum3, 'yyyy-mm-ddTHH:MM:SSZ');
    temp=[ge_quiver3(Long(n),Lat(n),0,u(n),v(n),w(n),...
                            'modelLinkStr',arrowStr,...
                            'altitudeMode','absolute',...
                            'arrowScale',5e1,...
                            'timeSpanStart', char(dstr), ...
                            'timeSpanStop', char(dstr2 ), ...
                            'name', 'quiver3 - red') ];
    k = findstr('longitude',temp);
    temp1 = [temp(1:k(1)+9) num2str(Long(n),10)  temp(k(2)-2:end)];
    l = findstr('latitude',temp1);
    temp2 = [temp1(1:l(1)+8) num2str(Lat(n),10)  temp1(l(2)-2:end)];
    
    s_red = [s_red temp2];
    
    
end


s_redf = ge_folder('red',s_red);

ge_output(kmlFileName,s_redf)
    
    