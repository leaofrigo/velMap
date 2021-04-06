function [VelExtrap] = FundoExtrapo(VelVec,metodo,ComecoCell,Depth,CellSize,NumCells,varargin)
if length(varargin)==2
    power = varargin{1};
    expsetting = varargin{2};
elseif length(varargin)==1
    power = varargin;
    expsetting.Z0 = 0.2; 
elseif isempty(varargin)==1
    power = 1/6;
    expsetting.Us = 0.1; 
end
ultimacellterm = ComecoCell+CellSize*NumCells;
highcell = abs(Depth-ultimacellterm);
xx = ultimacellterm:CellSize:Depth;
switch metodo
    case 'Constante'
        VelExtrap = ones(size(xx))*VelVec(end);
    case 'Exponencial'
        cells=linspace(0,highcell,NumCells);
        VelExtrap = expsetting.Us.*9.5.*(cells/expsetting.Z0).^power;        
    case 'No Slip'
        a = round(.2*length(VelVec));
        VelVect = VelVec(end-a:end);
        x = [ultimacellterm-a*CellSize:CellSize:ultimacellterm]';
        f=fit(x,VelVect,'a*exp(x)-a');
        VelExtrap = f(xx);
end
        
        

