function [VelU,VelL,VelR,temp1]= CalcularVelLVelR(ComecoCell,CellSize,NumOfCells,velN,velE,velU,Sumnum,th,VelU,VelL,VelR)
temp1 = -ComecoCell:-CellSize:-NumOfCells*CellSize-ComecoCell;
temp1 = temp1';
VelNTemp = velN(~isnan(velN));
VelETemp = velE(~isnan(velE));
VelUTemp = velU(~isnan(velU));
VelU(Sumnum(1):Sumnum(2)-1) = VelUTemp;
VelL(Sumnum(1):Sumnum(2)-1) = VelETemp*cos(th) - VelNTemp*sin(th);
VelR(Sumnum(1):Sumnum(2)-1) = VelETemp*sin(th) + VelNTemp*cos(th);   
end