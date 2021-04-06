clc; close all; clear all;
Directory = 'C:\Users\Roberta\Desktop\Ricardo\secoes\';
A = load ([Directory '20141107181441.mat']);

velocidades = A.WaterTrack.Velocity;
depth = - A.Summary.Depth;

%Calcular distancia viajada pelo barco
Dist = A.Summary.Track;
y = diff(Dist); % calcular diferenca entre pontos
y= sum(abs(y)); % somar a distancia absoluta entre pontos
largura=norm(y); % Achar a magnitude

vel_E = velocidades(:,1,:);
vel_N = velocidades(:,2,:);
vel_U = velocidades(:,3,:);
vel_D = velocidades(:,4,:);

TamanhoMatrix = size(velocidades);


for i = 1:TamanhoMatrix(3)-1
    velE(:,i) = vel_E(:,1,i);
    velN(:,i) = vel_N(:,1,i);
    velU(:,i) = vel_U(:,1,i);
    velD(:,i) = vel_D(:,1,i);
end


x = linspace(0, largura, length(depth));
ComecoCell = A.System.Cell_Start;
CellSize = A.System.Cell_Size;
NumOfCells  = A.Summary.Cells;

for k = 1:TamanhoMatrix(3)-1
%     if ceil(ComecoCell(k)*100)/100 == ComecoCell(k)
    temp = [-ComecoCell(k):-CellSize(k):-NumOfCells(k)*CellSize(k)-ComecoCell(k)]';
%     else
%         temp = [-ceil(ComecoCell(k)*100)/100:-CellSize(k):-NumOfCells(k)*CellSize(k)-ComecoCell(k)]';
%     end
    figure(1)
    hold on
    VelNTemp = velN(:,k);
    VelNTemp = VelNTemp(~isnan(VelNTemp));
    imagesc([x(k)],[ComecoCell(k),-temp(end)],VelNTemp);
    
    figure(2)
    hold on
    VelETemp = velE(:,k);
    VelETemp = VelETemp(~isnan(VelETemp));
    imagesc([x(k)],[ComecoCell(k),-temp(end)],VelETemp);
    
    figure(3)
    hold on
    VelUTemp = velU(:,k);
    VelUTemp = VelUTemp(~isnan(VelUTemp));
    imagesc([x(k)],[ComecoCell(k),-temp(end)],VelUTemp);
    
    figure(4)
    hold on
    VelDTemp = velD(:,k);
    VelDTemp = VelDTemp(~isnan(VelDTemp));
    imagesc([x(k)],[ComecoCell(k),-temp(end)],VelDTemp);
    
    figure(5)
    hold on
    VelMagTemp = sqrt(VelNTemp.^2+VelETemp.^2+VelUTemp.^2);
    imagesc([x(k)],[ComecoCell(k),-temp(end)],VelMagTemp);
%     quiver(x(k),-temp(1:end-1),-VelETemp,-VelNTemp)
    

end

for kk=1:5;
    figure(kk)
    colorbar
    caxis auto
%     caxis([-1 1]); %15-> [-3 1]
    t = colorbar('peer',gca); 
    set(get(t,'ylabel'),'String', 'Velocidade(m/s)');
    plot(x, -depth)
    xlim([0,largura])
    ylim([0,6])
    xlabel('Largura (m)')
    ylabel('Profundidade (m)')
    set(gca,'XDir','reverse')
    set(gca,'YDir','reverse')
    hold off
    
end
figure(1)
title(['Velocidade na direção Norte']);
figure(2)
title(['Velocidade na direção Leste']);
figure(3)
title(['Velocidade para cima']);
figure(4)
title(['Diferenca de velocidade (Erro de Velocidade)']);
figure(5)
title(['Magnitude da velocidade']);
% caxis([0 1.5]);
caxis auto;

% contourf(velU)
% colorbar; caxis([0 0.1180]); set(gca,'YDir','reverse') %15-> [-3 1]
% t = colorbar('peer',gca); 
% set(get(t,'ylabel'),'String', 'Velocidade(m/s)');


% fileID = fopen('arquivo2.dat');          
% A = textscan(fileID,'%f %f %f');                  
% fclose(fileID);                           
% cel = A{1};
% col = A{2};
% vel = A{3};




