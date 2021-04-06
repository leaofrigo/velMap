function distance = deg2distance(Lat,Long)
lat1=Lat(1);
lon1=Long(1);
lat2 = Lat(2);
lon2 = Long(2);
latrad1 = lat1*pi/180;
lonrad1 = lon1*pi/180;
latrad2 = lat2*pi/180;
lonrad2 = lon2*pi/180;

londif = abs(lonrad2-lonrad1);

raddis = acos(sin(latrad2)*sin(latrad1)+ ...
    cos(latrad2)*cos(latrad1)*cos(londif));
nautdis = raddis * 3437.74677;
stdiskm = nautdis * 1.852;
distance = stdiskm*1000;