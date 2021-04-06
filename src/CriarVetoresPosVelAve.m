function [Sumnum,Sumnum2,X,Y,AverageVelocityVector,AverageVelocityMag,cellDepths,zz,cellDepthsCentro,...
            VelETemp,VelNTemp,VelUTemp,VelDTemp,VelMagTemp,...
            Quiver3VectorX,Quiver3VectorY,Quiver3VectorZ,temp] =...
            CriarVetoresPosVelAve(ComecoCell,CellSize,NumOfCells,cellDepths,cellDepthsCentro,xx,...
            velN,velE,velU,velD,Sumnum,Sumnum2,X,Y,Quiver3VectorX,Quiver3VectorY,Quiver3VectorZ,Lat,Long)
% VelETemp,VelNTemp,VelUTemp,VelDTemp,VelMagTemp,
temp = -ComecoCell:-CellSize:-NumOfCells*CellSize-ComecoCell;
temp = temp';
cellDepths = [cellDepths;temp];
tempCentro = zeros(1,length(temp)-1);
for jj = 1:length(tempCentro);
     tempCentro(jj) = mean([temp(jj) temp(jj+1)]);
end
cellDepthsCentro = [cellDepthsCentro,tempCentro];
zz = abs(diff([temp(1),temp(end)]));
VelNTemp = velN(~isnan(velN));
VelETemp = velE(~isnan(velE));
VelUTemp = velU(~isnan(velU));
VelDTemp = velD(~isnan(velD));
VelMagTemp = sqrt(VelNTemp.^2+VelETemp.^2+VelUTemp.^2);
Sumnum = Sumnum+length(VelUTemp);
Sumnum2 = Sumnum2+length(temp);
X = [X xx*ones(1,length(VelUTemp))];
Y = [Y;temp(1:end-1)];
AverageVelocityVector = [mean(VelETemp) mean(VelNTemp) mean(VelUTemp)];
AverageVelocityMag = mean(VelMagTemp);
Quiver3VectorX = [Quiver3VectorX; ones(length(tempCentro),1)*Long];
Quiver3VectorY = [Quiver3VectorY; ones(length(tempCentro),1)*Lat];
Quiver3VectorZ = [Quiver3VectorZ; tempCentro(1:end)'];

end

