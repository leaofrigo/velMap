function [Area,depthnew,tracknew] = edges1(EdgeDist,depth,type,track,Edge01,handles)
numOfPoints = get(handles.Options2D.extrapnumpoints,'UserData');
% numOfPoints = 10;
if type ==1
    Area = 0;
    depthnew = depth;
    tracknew = depth;
    Disp('As aproximações laterais foram escolidas pelo usuario. Não foi possivel calcular extrapolação lateral')
else
    if Edge01 ==0; %Start Edge
        if type==1 %Vertical Bank
            Area = abs(depth(1)*EdgeDist);
            tracknew = [linspace(-EdgeDist,0,numOfPoints) track]; 
            depthnew = [ones(numOfPoints,1)*depth(1);depth];
        else %Sloped Bank,
            Area = abs(depth(1)*EdgeDist/2);
            tracknew = [linspace(track(1)-EdgeDist,track(1),numOfPoints) track];
            depthnew = [linspace(0,depth(1),numOfPoints)';depth];

        end
    else %End Edge
        if type ==1 %Vertical Bank
            Area = abs(depth(end)*EdgeDist);
            tracknew = [track linspace(track(end),track(end)+EdgeDist,numOfPoints)];
            depthnew = [depth ones(numOfPoints,1)*depth(end)];   
        else %Sloped Bank
            Area = abs(depth(end)*EdgeDist/2);
            tracknew = [track linspace(track(end),track(end)+EdgeDist,numOfPoints)]; 
            depthnew = [depth;linspace(depth(end),0,numOfPoints)'];

        end


    end
end
