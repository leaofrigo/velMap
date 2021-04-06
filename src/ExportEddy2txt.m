function output = ExportEddy2txt(Quiver3VectorXTot,Quiver3VectorYTot,...
    Quiver3VectorZTot,EDDYVECTOR,...
    YrR,XrR,FrR,handles)
Space = get(handles.Options3D.NumofPointsBat,'UserData');
data = get(handles.Options3D.SimulOpt,'UserData');
xmin = min(Quiver3VectorXTot(:)); 
ymin = min(Quiver3VectorYTot(:)); 
zmin = min(Quiver3VectorZTot(:));

xmax = max(Quiver3VectorXTot(:)); 
ymax = max(Quiver3VectorYTot(:)); 
zmax = max(Quiver3VectorZTot(:));

x1 = linspace(xmin,xmax,Space);
y1 = linspace(ymin,ymax,Space);
z1 = linspace(zmin,zmax,Space);
[X1,Y1,Z1] = meshgrid(x1,y1,z1);
ss= size(X1);
VEL3D = griddata(Quiver3VectorXTot,Quiver3VectorYTot,Quiver3VectorZTot,EDDYVECTOR,X1,Y1,Z1);
Vel = EliminateVSmallerThanBat(X1,Y1,Z1,VEL3D,XrR,YrR,FrR);
NumOfLayers = handles.EddyExpLayers;
LayerDistance = handles.EddyDistLayers;
LayerLevels = max(Z1(:)):-LayerDistance:-LayerDistance*NumOfLayers+max(Z1(:))+.1;
LayerXT = nan(NumOfLayers,ss(1),ss(2));
LayerYT = LayerXT;
EddyT = LayerXT;
LayerZT = ones(NumOfLayers,ss(1),ss(2));
Dir = [handles.Directory 'EddyViscosityExport\'];
A = exist(Dir,'dir');
if A
else
    mkdir(Dir);
end
for i = 1 : length(LayerLevels)
    if i ~= length(LayerLevels)
        tf = and(Z1(1,1,:)<LayerLevels(i),Z1(1,1,:)>LayerLevels(i+1));
        tf = tf(:);
        LayerX = X1(:,:,tf);
        LayerXT(i,:,:) = nanmean(LayerX,3);
        LayerY = Y1(:,:,tf);
        LayerYT(i,:,:) = nanmean(LayerY,3);
        Eddy = Vel(:,:,tf);
        EddyT(i,:,:) = nanmean(Eddy,3); 
        LayerZT(i,:,:) = LayerZT(i,:,:)*mean(LayerLevels(i:i+1));
        
        
%         tf = and(Z1<LayerLevels(i),Z1>LayerLevels(i+1));
%         LayerX = X1(tf);
%         LayerY = Y1(tf);
%         Eddy = Vel(tf);
%         LayerZ = ones(size(Eddy))*mean(LayerLevels(i:i+1));
%         

    else
        tf = Z1(1,1,:)<LayerLevels(i);
        tf = tf(:);
        LayerX = X1(:,:,tf);
        LayerXT(i,:,:) = nanmean(LayerX,3);
        LayerY = Y1(:,:,tf);
        LayerYT(i,:,:) = nanmean(LayerY,3);
        Eddy = Vel(:,:,tf);
        EddyT(i,:,:) = nanmean(Eddy,3); 
        LayerZT(i,:,:) = LayerZT(i,:,:)*(-LayerDistance*NumOfLayers+max(Z1(:)));
%         tf = Z1<LayerLevels(i);
%         LayerX = X1(tf);
%         LayerY = Y1(tf);
%         Eddy = Vel(tf);
%         LayerZ= ones(size(Eddy))*(-LayerDistance*NumOfLayers+max(Z1(:)));
    end
    File = fopen([Dir 'Layer_' num2str(i) '.txt'],'w+');
    formatSpec = '%3.8f; %3.8f; %4.1f; %1.4e;\n';
    LayerXR = reshape(LayerXT,1,[]);
    LayerYR = reshape(LayerYT,1,[]);
    LayerZR = reshape(LayerZT,1,[]);
    EddyR = reshape(EddyT,1,[]);
    TF = ~isnan(EddyR);
%     size(LayerXR)
%     size(LayerYR)
%     size(LayerZR)
%     size(EddyR)
%     size(TF)
%     LayerXR(TF);
%     LayerYR(TF);
%     LayerZR(TF);
%     EddyR(TF);
    fprintf(File,formatSpec,[LayerXR(TF);LayerYR(TF); LayerZR(TF); EddyR(TF)]);
    
    assignin('base','LayerXT',LayerXT)
    assignin('base','LayerYT',LayerYT)
    assignin('base','LayerZT',LayerZT)
    assignin('base','EddyT',EddyT)
    fclose(File);
end
output.LayerX = LayerXT;
output.LayerY = LayerYT;
output.LayerZ = LayerZT;
output.Value = EddyT;
% for j =1:NumOfLayers
%     File = fopen([Dir 'Layer_' num2str(j) '.txt'],'w+');
% %     FileY = fopen([Dir 'LayerLatitude_' num2str(j) '.txt'],'w+');
% %     FileZ = fopen([Dir 'LayerDepth_' num2str(j) '.txt'],'w+');
% %     FileEddy = fopen([Dir 'LayerEddy_' num2str(j) '.txt'],'w+');
%     for k=1:length(LayerXT)
%         fprintf(FileX,'%3.8f; ',LayerXT(j,k,:));
%         fprintf(FileX,'\n');
% %         fprintf(FileY,'%3.8f; ',LayerYT(j,k,:));
% %         fprintf(FileY,'\n');
% %         fprintf(FileZ,'%4.1f; ',LayerZT(j,k,:));
% %         fprintf(FileZ,'\n');
% %         fprintf(FileEddy,'%10.0f; ',EddyT(j,k,:));
% %         fprintf(FileEddy,'\n');
%     end
%     fclose(FileX);
%     fclose(FileY);
%     fclose(FileZ);
%     fclose(FileEddy);
% end
    
    
