function writeExcelReport(excel,path)
titles = {'Section:';'Total Q (m^3/s)';'Mean Water Speed (m/s)';'Distance Traveled (m)';...
    'Mag. Distance Traveled (m)';'Max depth Cell Start (m)';...
    'Perimeter (m)';'Area Total (m^2)';'Hydraulic Radius (m)';'Average Height (m)';'Froude number (Fr)';...
    'Area Measured/Area Total'};
AveHeight = zeros(length(excel));
FrNum = zeros(length(excel));
xlswrite(path,titles,'All sessions','A1')
for i = 1:length(excel)
    AveHeight(i) = excel(i).AreaTot/excel(i).diffxx1xxend;
    FrNum(i) =  excel(i).Ave_speed/sqrt(9.80665*AveHeight(i));
    sheet=excel(i).section;
    values = [[sheet{1} '.mat'];{excel(i).Total_Q};{excel(i).Ave_speed};{excel(i).totaldist};...
        {excel(i).diffxx1xxend};{excel(i).maxdepthstart};{excel(i).perimeter};{excel(i).AreaTot};...
        {excel(i).hidrrad};{AveHeight(i)};{FrNum(i)};{excel(i).AreaMeas/excel(i).AreaTot}];
    All = [titles values];

    xlswrite(path,All,[sheet{1} '.mat'],'A1')
    column = num2letter(i+1);
    xlswrite(path,values,'All sessions',[column '1'])
end

