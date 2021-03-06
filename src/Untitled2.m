k = kml('Contourf');
f = k.createFolder('kml.contourf');
f.contourf(x,y,z,'description','Eddy Viscosity at surface','numberOfLevels',10);


k1 = kml('Contour');
f1 = k1.createFolder('kml.contour');
f1.contour(x,y,z,'description','Eddy Viscosity at surface','numberOfLevels',10,'altitude',2000,'altitudeMode','relativeToGround')



k2 = kml('Transfer');

fh = figure;
ax = gca;

[~,h]= contourf(ax,x,y,z,10);
set(h,'LineColor','none')
set(h,'LineStyle','none')
% set(h,'ShowText','on')
set(h,'HandleVisibility','off')
xlabel('longitude');
ylabel('latitude');
snapnow;
k2.transfer(ax,'keepAxis',false,'transparentBG',true,'altitudeMode','absolute','altitude',1000);
k.run;
k1.run;
k2.run
close(fh)