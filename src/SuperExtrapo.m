function [VelExtrap] = SuperExtrapo(VelVec,metodo,ComecoCell,CellSize,varargin)
if length(varargin)==1
    power = varargin;
elseif isempty(varargin)==1
    power = 1/6;
else
    error('No maximo cinco inputs. (Vetor de Vel, Metodo, Comeco da Celula,Tamanho da Celula,Valor Exponencial)')
end
xx = 0:CellSize:ComecoCell;
switch metodo
    case 'Constante'
        VelExtrap = ones(size(xx))*VelVec(1);
    case 'Exponencial'
        cells=linspace(0,highcell,CelNum);
        VelExtrap = expsetting.Us.*9.5.*(cells/expsetting.Z0).^power;

        
    case 'Tres Pontos'
        x =[ComecoCell:CellSize:ComecoCell+CellSize*2];
        y =VelVec(1:3);
        f = fit(x,y,'poly1');
        VelExtrap = f(xx);
end