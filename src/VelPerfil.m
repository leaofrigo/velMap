function [AverageVel, CellDepth] = VelPerfil(CellSize,CellStart,Sumnum,Vel,NumOfCellsExit,xx)
NumOfCellPCol = diff(Sumnum)';
IndDepth = CellStart+CellSize.*NumOfCellPCol;
IndLength = diff(xx)';
MaxDepth = max(IndDepth);
MinDepth = min(CellStart);
Steps = linspace(MinDepth,MaxDepth,NumOfCellsExit+1);
AverageVel = zeros(NumOfCellsExit,1);
TF = false([NumOfCellsExit length(Vel)]);
AllCellDepth = zeros(size(Vel));
AllCellSize = zeros(size(Vel));
AllCellLength = zeros(size(Vel));
CellDepth = zeros(NumOfCellsExit,1);
for m = 1:length(IndLength)
    as = CellStart(m)+.5*CellSize(m):CellSize(m):CellSize(m)*NumOfCellPCol(m)+CellStart(m)-.5*CellSize(m)';
    AllCellSize(Sumnum(m)+1:Sumnum(m+1)) = ones(size(as))*CellSize(m);
    AllCellLength(Sumnum(m)+1:Sumnum(m+1)) = ones(size(as))*IndLength(m);
    AllCellDepth(Sumnum(m)+1:Sumnum(m+1)) = as;   
end
for n =1:NumOfCellsExit    
    TF(n,:) = and(AllCellDepth>=Steps(n),AllCellDepth<Steps(n+1));
    CellDepth(n) = mean([Steps(n) Steps(n+1)]);
    AverageVel(n) = nansum(Vel(TF(n,:)).*AllCellSize(TF(n,:)).*AllCellLength(TF(n,:)))./nansum(AllCellSize(TF(n,:)).*AllCellLength(TF(n,:)));
end







