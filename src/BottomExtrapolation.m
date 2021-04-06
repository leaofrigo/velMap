function [VelMagExt,VelEExt,VelNExt,VelUExt,VelDExt] = BottomExtrapolation(Type,Profile,Power,NumOfCells,VelMag,VelE,VelN,VelU,VelD,CellHeight,depth)
if Type == 0 %Constant Fit
    if Profile == 0 %use entire profile
        VelMagExt = mean(VelMag);
        VelEExt = mean(VelE);
        VelNExt = mean(VelN);
        VelUExt = mean(VelU);
        VelDExt = mean(VelD);
    else %user Selected Number of Cells
        VelMagExt = mean(VelMag(end-NumOfCells+1:end))
        VelEExt = mean(VelE(end-NumOfCells+1:end));
        VelNExt = mean(VelN(end-NumOfCells+1:end))
        VelUExt = mean(VelU(end-NumOfCells+1:end))
        VelDExt = mean(VelD(end-NumOfCells+1:end))      
    end   
else %Power Fit
    numofextcells = 5;
    if Profile == 0 %use entire profile
        fVelMagExt=fit(VelMag,CellHeight(2:end),['b+ a*exp(x*' mat2str(Power) ')']);
        VelMagExt = fVelMagExt(linspace(CellHeight(end),depth,numofextcells));
        fVelEExt=fit(VelE,CellHeight(2:end),['b+ a*exp(x*' mat2str(Power) ')']);
        VelEExt = fVelEExt(linspace(CellHeight(end),depth,numofextcells));
        fVelNExt=fit(VelN,CellHeight(2:end),['b+ a*exp(x*' mat2str(Power) ')']);
        VelNExt = fVelNExt(linspace(CellHeight(end),depth,numofextcells));
        fVelUExt=fit(VelU,CellHeight(2:end),['b+ a*exp(x*' mat2str(Power) ')']);
        VelUExt = fVelUExt(linspace(CellHeight(end),depth,numofextcells));
        fVelDExt=fit(VelD,CellHeight(2:end),['b+ a*exp(x*' mat2str(Power) ')']);
        VelDExt = fVelDExt(linspace(CellHeight(end),depth,numofextcells));      
    else %user Selected Number of Cells
        fVelMagExt=fit(VelMag(end-NumOfCells+1:end),CellHeight(end-NumOfCells+1:end),['b+ a*exp(x*' mat2str(Power) ')']);
        VelMagExt = fVelMagExt(linspace(CellHeight(1),0,numofextcells));
        fVelEExt=fit(VelE(end-NumOfCells+1:end),CellHeight(end-NumOfCells+1:end),['b+ a*exp(x*' mat2str(Power) ')']);
        VelEExt = fVelEExt(linspace(CellHeight(1),0,numofextcells));
        fVelNExt=fit(VelN(end-NumOfCells+1:end),CellHeight(end-NumOfCells+1:end),['b+ a*exp(x*' mat2str(Power) ')']);
        VelNExt = fVelNExt(linspace(CellHeight(1),0,numofextcells));
        fVelUExt=fit(VelU(end-NumOfCells+1:end),CellHeight(end-NumOfCells+1:end),['b+ a*exp(x*' mat2str(Power) ')']);
        VelUExt = fVelUExt(linspace(CellHeight(1),0,numofextcells));
        fVelDExt=fit(VelD(end-NumOfCells+1:end),CellHeight(end-NumOfCells+1:end),['b+ a*exp(x*' mat2str(Power) ')']);
        VelDExt = fVelDExt(linspace(CellHeight(1),0,numofextcells));    
    end    
end
end