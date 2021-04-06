function LayersKML(Input,numColors)

X = Input.LayerX;
Y = Input.LayerY;
Z = Input.LayerZ;
V = Input.Value;
Alt = Input.Alt;

S = size(V);
V = abs(V);
fh = figure;
ax = gca;
% [~, UpperLim,~] = limites(V);
for i = 1:S(1)+1
    k(i) = kml(['Contourf Layer ' num2str(i)]);
    
    Xtemp = reshape(X(i,:,:),S(2),S(3));
    Ytemp = reshape(Y(i,:,:),S(2),S(3));
    Ztemp = reshape(Z(i,:,:),1,[]);
    Vtemp = reshape(V(i,:,:),S(2),S(3));
    
    [~,h]= contourf(ax,Xtemp,Ytemp,Vtemp,numColors);

    caxis([0,.03])
    set(h,'LineColor','none')
    set(h,'LineStyle','none')
    set(h,'ShowText','on')
%     set(h,'HandleVisibility','off')
%     xlabel('longitude');
%     ylabel('latitude');
    snapnow;
    k(i).transfer(ax,'keepAxis',false,'transparentBG',true,'altitudeMode',...
        'absolute','altitude',Alt+300*(S(1)-i));
end
for i = 1:S(1)
    k(i).run
end
close(fh)
    