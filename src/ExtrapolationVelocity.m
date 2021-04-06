function [VelocityExtrapolated,Y] = ExtrapolationVelocity(ExtrapStruct)
VelVec = ExtrapStruct.VelVec;
s=size(VelVec);
s=s(2);
s1=s(1);
metodo = ExtrapStruct.metodo;
TopOrBottom = ExtrapStruct.type;
NumCellsExtrap = ExtrapStruct.NumCellsExtrap;
Depth = ExtrapStruct.Depth;
CellStart = ExtrapStruct.CellStart;

switch TopOrBottom
    case 'Top'
        Y = linspace(abs(-Depth-CellStart),-Depth,NumCellsExtrap);
        VelVecCon = VelVec(1);
    case 'Bottom'
        NumCellsInColumn = ExtrapStruct.NumCellsInColumn;
        CellSize = ExtrapStruct.CellSize;
        ultimacellterm = CellStart(1)+CellSize*NumCellsInColumn;
        Y = linspace(0,abs(-Depth-ultimacellterm),NumCellsExtrap);
        VelVecCon = VelVec(end);
end
if metodo == 1
    metodo = 'Exponencial';
elseif metodo ==0
    metodo = 'Constante';
end
Y=Y';
VelocityExtrapolated= zeros(s,length(Y));

switch metodo
    case 'Constante'
        VelocityExtrapolated = ones(size(Y))*VelVecCon;
    case 'Exponencial'
        power = ExtrapStruct.power;
        Z0 = ExtrapStruct.expsetting.Z0;
        Ustar = ExtrapStruct.expsetting.Ustar;
        VelocityExtrapolated = Ustar.*9.5.*(Y/Z0).^power;        
    case 'No Slip'
        a = round(.2*s1);
        VelVectexp = VelVec(end-a:end);
        x = ultimacellterm-a*CellSize:CellSize:ultimacellterm;
        f=fit(x,VelVectexp,'a*exp(x)-a');
        VelocityExtrapolated = f(Y);
    case 'Tres Pontos'
        x =CellStart:CellSize:CellStart+CellSize*2;
        y =VelVec(1:3);
        f = fit(x,y,'poly1');
        VelocityExtrapolated = f(Y);
end
