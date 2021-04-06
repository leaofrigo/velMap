clear all;close all; clc

latlim = [38.895003  38.900104];
lonlim = [-77.041641 -77.030537];

numberOfAttempts = 5;
attempt = 0;
info = [];
serverURL = 'http://raster.nationalmap.gov/arcgis/services/Orthoimagery/USGS_EROS_Ortho/ImageServer/WMSServer?';
while(isempty(info))
    try
        info = wmsinfo(serverURL);
        orthoLayer = info.Layer(1);
    catch e 
        
        attempt = attempt + 1;
        if attempt > numberOfAttempts
            throw(e);
        else
            fprintf('Attempting to connect to server:\n"%s"\n', serverURL)
        end        
    end
end

imageLength = 1024;
[A, R] = wmsread(orthoLayer, 'Latlim', latlim, 'Lonlim', lonlim, ...
    'ImageHeight', imageLength, 'ImageWidth', imageLength);

figure
axesm('utm', 'Zone', utmzone(latlim, lonlim), ...
    'MapLatlimit', latlim, 'MapLonlimit', lonlim, ...
    'Geoid', wgs84Ellipsoid)
geoshow(A,R)
axis off
% title({'San Francisco','Northern Section of Golden Gate Bridge'})