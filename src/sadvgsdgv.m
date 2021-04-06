figure('units','normalized');
%plot first contour in axes1
Z1 = peaks;
ax(1)=axes;
[C1,h1] = contour(Z1);
hc(1)=colorbar;
ax(2)=axes;
Z2=rot90(interp2(10*Z1,2));
[C2,h2] = contour(Z2);
%move x- and y-axis of axes2
set(ax(2), 'XAxisLocation','top',...
             'YAxisLocation','right',...
             'Color','none');      
hc(2)=colorbar;
%reduce axes position to make room for the second colorbar
pos=get(ax(1),'Position');
set(ax,'Position',[pos(1) pos(2) 0.9*pos(3) pos(4)]);
cpos=get(hc(1),'position');
%move the second colorbar to the right
set(hc(2),'position',[cpos(1)+0.05 cpos(2:4)]);
