function [TopVelMagExt,TopVelEExt,TopVelNExt,TopVelUExt,TopVelDExt,...
    BottomVelMagExt,BottomVelEExt,BottomVelNExt,BottomVelUExt,BottomVelDExt,numcel]...
    = TopBottomExtrapolation(Type,Profile,Power,NumOfCells,VelMag,VelE,VelN,VelU,VelD,CellHeight,depth)
numcel = 5;
if Type == 1 %Constant Fit 'Default = 0'
    if Profile == 0 %use entire profile
        TopVelMagExt = ones(numcel,1)*mean(VelMag);
        BottomVelMagExt = TopVelMagExt;
        TopVelEExt = ones(numcel,1)*mean(VelE);
        BottomVelEExt = TopVelEExt;
        TopVelNExt = ones(numcel,1)*mean(VelN);
        BottomVelNExt = TopVelNExt;
        TopVelUExt = ones(numcel,1)*mean(VelU);
        BottomVelUExt = TopVelUExt;
        TopVelDExt = ones(numcel,1)*mean(VelD);
        BottomVelDExt = TopVelDExt;
    else %user Selected Number of Cells
        TopVelMagExt = ones(numcel,1)*mean(VelMag(end-NumOfCells+1:end));
        BottomVelMagExt = TopVelMagExt;
        TopVelEExt = ones(numcel,1)*mean(VelE(end-NumOfCells+1:end));
        BottomVelEExt = TopVelEExt;
        
        TopVelNExt = ones(numcel,1)*mean(VelN(end-NumOfCells+1:end));
        BottomVelNExt = TopVelNExt;
        TopVelUExt = ones(numcel,1)*mean(VelU(end-NumOfCells+1:end));
        BottomVelUExt = TopVelUExt;
        TopVelDExt = ones(numcel,1)*mean(VelD(end-NumOfCells+1:end)) ;
        BottomVelDExt = TopVelDExt;
    end   
else %Power Fit
    FitLine = ['b+a*exp(x*' mat2str(Power) ')'];
    if Profile == 0 %use entire profile
        fVelMagExt=fit(VelMag,CellHeight(2:end),FitLine);
        TopVelMagExt = fVelMagExt(linspace(CellHeight(1),0,numcel));
        BottomVelMagExt = fVelMagExt(linspace(CellHeight(end),depth,numcel));
        fVelEExt=fit(VelE,CellHeight(2:end),FitLine);
        TopVelEExt = fVelEExt(linspace(CellHeight(1),0,numcel));
        BottomVelEExt = fVelEExt(linspace(CellHeight(end),depth,numcel));
        fVelNExt=fit(VelN,CellHeight(2:end),FitLine);
        TopVelNExt = fVelNExt(linspace(CellHeight(1),0,numcel));
        BottomVelNExt = fVelNExt(linspace(CellHeight(end),depth,numcel));
        fVelUExt=fit(VelU,CellHeight(2:end),FitLine);
        TopVelUExt = fVelUExt(linspace(CellHeight(1),0,numcel));
        BottomVelUExt = fVelUExt(linspace(CellHeight(end),depth,numcel));
        fVelDExt=fit(VelD,CellHeight(2:end),FitLine);
        TopVelDExt = fVelDExt(linspace(CellHeight(1),0,numcel));
        BottomVelDExt = fVelDExt(linspace(CellHeight(end),depth,numcel));      
    else %user Selected Number of Cells
        fVelMagExt=fit(VelMag(end-NumOfCells+1:end),CellHeight(end-NumOfCells+1:end),FitLine);
        TopVelMagExt = fVelMagExt(linspace(CellHeight(1),0,numcel));
        BottomVelMagExt = fVelMagExt(linspace(CellHeight(1),0,numcel));
        fVelEExt=fit(VelE(end-NumOfCells+1:end),CellHeight(end-NumOfCells+1:end),FitLine);
        TopVelEExt = fVelEExt(linspace(CellHeight(1),0,numcel));
        BottomVelEExt = fVelEExt(linspace(CellHeight(1),0,numcel));
        fVelNExt=fit(VelN(end-NumOfCells+1:end),CellHeight(end-NumOfCells+1:end),FitLine);
        TopVelNExt = fVelNExt(linspace(CellHeight(1),0,numcel));
        BottomVelNExt = fVelNExt(linspace(CellHeight(1),0,numcel));
        fVelUExt=fit(VelU(end-NumOfCells+1:end),CellHeight(end-NumOfCells+1:end),FitLine);
        TopVelUExt = fVelUExt(linspace(CellHeight(1),0,numcel));
        BottomVelUExt = fVelUExt(linspace(CellHeight(1),0,numcel));
        fVelDExt=fit(VelD(end-NumOfCells+1:end),CellHeight(end-NumOfCells+1:end),FitLine);
        TopVelDExt = fVelDExt(linspace(CellHeight(1),0,numcel)); 
        BottomVelDExt = fVelDExt(linspace(CellHeight(1),0,numcel));    
    end    
end
end