function [VelMagExt,VelEExt,VelNExt,VelUExt,VelDExt] = TopExtrapolation(Type,Profile,Power,NumOfCells,VelMag,VelE,VelN,VelU,VelD,CellHeight)
if Type == 0 %Constant Fit
    if Profile == 0 %use entire profile
        VelMagExt = mean(VelMag);
        VelEExt = mean(VelE);
        VelNExt = mean(VelN);
        VelUExt = mean(VelU);
        VelDExt = mean(VelD);
    else %user Selected Number of Cells
        VelMagExt = mean(VelMag(1:NumOfCells));
        VelEExt = mean(VelE(1:NumOfCells));
        VelNExt = mean(VelN(1:NumOfCells));
        VelUExt = mean(VelU(1:NumOfCells));
        VelDExt = mean(VelD(1:NumOfCells));       
    end
else %Power Fit
    numofextcells = 5;
    if Profile == 0 %use entire profile
        fVelMagExt=fit(VelMag,CellHeight(1:end-1),['a*exp(x*' mat2str(Power) ')']);
        VelMagExt = fVelMagExt(linspace(CellHeight(1),0,numofextcells));
        fVelEExt=fit(VelE,CellHeight(1:end-1),['a*exp(x*' mat2str(Power) ')']);
        VelEExt = fVelEExt(linspace(CellHeight(1),0,numofextcells));
        fVelNExt=fit(VelN,CellHeight(1:end-1),['a*exp(x*' mat2str(Power) ')']);
        VelNExt = fVelNExt(linspace(CellHeight(1),0,numofextcells));
        fVelUExt=fit(VelU,CellHeight(1:end-1),['a*exp(x*' mat2str(Power) ')']);
        VelUExt = fVelUExt(linspace(CellHeight(1),0,numofextcells));
        fVelDExt=fit(VelD,CellHeight(1:end-1),['a*exp(x*' mat2str(Power) ')']);
        VelDExt = fVelDExt(linspace(CellHeight(1),0,numofextcells));       
    else %user Selected Number of Cells
        fVelMagExt=fit(VelMag(1:NumOfCells),CellHeight(1:NumOfCells),['a*exp(x*' mat2str(Power) ')']);
        VelMagExt = fVelMagExt(linspace(CellHeight(1),0,numofextcells));
        fVelEExt=fit(VelE(1:NumOfCells),CellHeight(1:NumOfCells),['a*exp(x*' mat2str(Power) ')']);
        VelEExt = fVelEExt(linspace(CellHeight(1),0,numofextcells));
        fVelNExt=fit(VelN(1:NumOfCells),CellHeight(1:NumOfCells),['a*exp(x*' mat2str(Power) ')']);
        VelNExt = fVelNExt(linspace(CellHeight(1),0,numofextcells));
        fVelUExt=fit(VelU(1:NumOfCells),CellHeight(1:NumOfCells),['a*exp(x*' mat2str(Power) ')']);
        VelUExt = fVelUExt(linspace(CellHeight(1),0,numofextcells));
        fVelDExt=fit(VelD(1:NumOfCells),CellHeight(1:NumOfCells),['a*exp(x*' mat2str(Power) ')']);
        VelDExt = fVelDExt(linspace(CellHeight(1),0,numofextcells));      
    end    
end
end