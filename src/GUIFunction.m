function varargout = GUIFunction(varargin)
% GUIFUNCTION MATLAB code for GUIFunction.fig
%      GUIFUNCTION, by itself, creates a new GUIFUNCTION or raises the existing
%      singleton*.
%
%      H = GUIFUNCTION returns the handle to a new GUIFUNCTION or the handle to
%      the existing singleton*.
%
%      GUIFUNCTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIFUNCTION.M with the given input arguments.
%
%      GUIFUNCTION('Property','Value',...) creates a new GUIFUNCTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUIFunction_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUIFunction_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUIFunction

% Last Modified by GUIDE v2.5 18-Jun-2015 17:27:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUIFunction_OpeningFcn, ...
                   'gui_OutputFcn',  @GUIFunction_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GUIFunction is made visible.
function GUIFunction_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUIFunction (see VARARGIN)

% Choose default command line output for GUIFunction
p = mfilename('fullpath');
UserSettings = load([p '\UserSettings.mat']);
UserSettings = UserSettings.UserSettings;
handles.UserSettings = UserSettings;
handles.run2d = 0;
handles.run3d = 0;
addpath([p(1:end-11) 'googleearth'])

handles.Options2D.Smoo2dOption = handles.SmoothMenu;
handles.Options2D.quiv2dOption = handles.QuiverMenu;
handles.Options2D.shpfileexp = handles.ShpFileMenu;
handles.Options2D.ExpShp2d = handles.ShapeType;
handles.Options2D.plotpostion = handles.PlotPositionMenu;
handles.Options2D.extrapnumpoints = handles.SideExtrapMenu;
handles.Options2D.VortSide =handles.SideViewMenu;
handles.Options2D.Z0Ustar = handles.ChenMenu;
handles.Options2D.NumCellsExtrap = handles.VelExtrapMenu;
handles.Options2D.KML2doption = handles.KmlMenu;

set(handles.Options2D.Smoo2dOption,'UserData',UserSettings.Options2D.Smoo2dValue);
set(handles.Options2D.quiv2dOption,'UserData',UserSettings.Options2D.quiver2dValue);
set(handles.Options2D.plotpostion,'UserData',UserSettings.Options2D.plotpostionValue);
set(handles.Options2D.ExpShp2d,'UserData',UserSettings.Options2D.ExpShp2dValue);
set(handles.Options2D.shpfileexp,'UserData',UserSettings.Options2D.shpfileexpValue);
set(handles.Options2D.extrapnumpoints,'UserData',UserSettings.Options2D.extrapnumpointsValue);
set(handles.Options2D.VortSide,'UserData',UserSettings.Options2D.VortSideValue);
set(handles.Options2D.Z0Ustar,'UserData',UserSettings.Options2D.Z0Ustar1Value);
set(handles.Options2D.NumCellsExtrap,'UserData',UserSettings.Options2D.NumCellsExtrapValue);
set(handles.Options2D.KML2doption,'UserData',UserSettings.Options2D.KML2dValue);


handles.Options3D.Smoo3dOption = handles.Smooth2DMenu3D;
handles.Options3D.quiv3dOption = handles.Quiver2DMenu3D;
handles.Options3D.ExpShp3d = handles.TypeShpMenu3D;
handles.Options3D.KML3doption = handles.KmlMenu3D;
handles.Options3D.quiver3doption = handles.Quiver3DMenu3D;
handles.Options3D.NumofPointsBat = handles.BathPointsMenu3D;
handles.Options3D.smallgraphlimts = handles.MiniPLotMenu3D;
handles.Options3D.plotpostion3d = handles.PlotPositionMenu3D;
handles.Options3D.shpfileexp3d = handles.FileShpMenu3D;
handles.Options3D.NumCellsExtrap = handles.VelExtrapMenu3D;
handles.Options3D.Z0Ustar = handles.ChenMenu3D;
handles.Options3D.SimulOpt = handles.AnimationMenu3D;
handles.Options3D.VisTopQuiv = handles.TopViewQuivMenu3D;
handles.Options3D.GridTopQuiv = handles.TopViewGridMenu3D;
handles.Options3D.ContourLines = handles.NumColorBathMenu3D;
handles.Options3D.ContourMap = handles.ColorMapMenu3D;


set(handles.Options3D.Smoo3dOption,'UserData',UserSettings.Options3D.Smoo2dValue);
set(handles.Options3D.quiv3dOption,'UserData',UserSettings.Options3D.quiver2dValue);
set(handles.Options3D.plotpostion3d,'UserData',UserSettings.Options3D.plotpostionValue);
set(handles.Options3D.ExpShp3d,'UserData',UserSettings.Options3D.ExpShp2dValue);
set(handles.Options3D.shpfileexp3d,'UserData',UserSettings.Options3D.shpfileexpValue);
set(handles.Options3D.Z0Ustar,'UserData',UserSettings.Options3D.Z0Ustar1Value);
set(handles.Options3D.NumCellsExtrap,'UserData',UserSettings.Options3D.NumCellsExtrapValue);
set(handles.Options3D.KML3doption,'UserData',UserSettings.Options3D.KML3dValue);
set(handles.Options3D.quiver3doption,'UserData',UserSettings.Options3D.extrapnumpointsValue);
set(handles.Options3D.smallgraphlimts,'UserData',UserSettings.Options3D.smallgraphlimts);
set(handles.Options3D.SimulOpt,'UserData',UserSettings.Options3D.SimulOptValue);
set(handles.Options3D.VisTopQuiv,'UserData',UserSettings.Options3D.VisTopQuivValue);
set(handles.Options3D.GridTopQuiv,'UserData',UserSettings.Options3D.GridTopQuivValue);
set(handles.Options3D.ContourLines,'UserData',UserSettings.Options3D.ContourLinesValue);
set(handles.Options3D.ContourMap,'UserData',UserSettings.Options3D.ContourMapValue);
set(handles.Options3D.quiver3doption,'UserData',UserSettings.Options3D.quiver3dValue);
set(handles.Options3D.NumofPointsBat,'UserData',UserSettings.Options3D.NumofPointsBatValue);
% 
% 
% 



handles.output = hObject;
Directory = uigetdir(pwd);

% Directory = uigetdir('C:\Users\Roberta\Desktop\Ricardo\detalhamento do tauri adcp\seçoes pra rodar\Novas\');
handles.Directory = [Directory '\'];
OpenFiles = DirectoryOpenFiles(Directory);
handles.OpenFiles = OpenFiles;
set(handles.PopSecoes,'String',OpenFiles)
GraficosPPlotar = [{'North'} {'East'} {'Upstream'} {'Magnitude'} {'Error'} {'Transverse'} {'Longitudinal'}];
handles.GraficosPPlotar = GraficosPPlotar;
set(handles.Graficos,'String',GraficosPPlotar)
ColumnFormat={'char','logical'};
ColumnName = {'Sessions','Process'};
data1=[];
for i = 1:length(OpenFiles)
    data1 = [data1;{OpenFiles{i}(1:end-4) true}];
end
set(handles.Secoes3D,'data',data1,'ColumnName',ColumnName,'ColumnFormat',ColumnFormat,...
    'ColumnEditable', [false true],'RowName',[],'ColumnWidth',{110, 75});
% Inicializar Variaveis
handles.plotpreview = 0;
handles.plotpreview2d = 0;
handles.GraficosV = 'North';
if length(OpenFiles)<1
    handles.PopSecoesV = '';
else
    handles.PopSecoesV = OpenFiles{1};
end
handles.SmoothV = 0;
handles.Kml2DV = 0;
handles.BatimetriaV = 0;
handles.GridV = 0;
handles.Kml3dV = 0;
handles.quiver3dV = 0;
handles.simulacaoV = 0;
handles.tabledata = data1;
handles.SuperiorV = 0;
handles.Excel2DV = 0;
handles.Excel3DV= 0;
handles.PlotSaveAllV = 0;
handles.ShapeV = 0;
handles.Shape3DV = 0;
handles.PlotV =1;
handles.quiver2dV =0;
handles.quiver2d3dV = 0;
handles.Smooth2d3dV = 0;
handles.StreamLinesV = 0;
handles.topVorticity = 0;
handles.SideVort.V = 0;
handles.exp2grid2D.V = 0;
DisplayDir_Callback(handles.DisplayDir, eventdata, handles)
% Update handles structure
guidata(hObject,handles)
% ResetDefault3d = uimenu(menu1,'Label','ResetDefault Options3D','Callback',{@Reset3D,handles});
% ResetDefault = uimenu(menu,'Label','ResetDefault Options2D','Callback',...
%     {@Reset2D,handles});
% UIWAIT makes GUIFunction wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUIFunction_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in Graficos.
function Graficos_Callback(hObject, eventdata, handles)
% hObject    handle to Graficos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = get(hObject,'Value');
GraficosPPlotar = handles.GraficosPPlotar;
handles.GraficosV = GraficosPPlotar{contents};
guidata(hObject,handles)
% Hints: contents = cellstr(get(hObject,'String')) returns Graficos contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Graficos


% --- Executes during object creation, after setting all properties.
function Graficos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Graficos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in PopSecoes.
function PopSecoes_Callback(hObject, eventdata, handles)
% hObject    handle to PopSecoes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = get(hObject,'Value');
OpenFiles = handles.OpenFiles;
handles.PopSecoesV = OpenFiles{contents};
guidata(hObject,handles)
% Hints: contents = cellstr(get(hObject,'String')) returns PopSecoes contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PopSecoes



% --- Executes during object creation, after setting all properties.
function PopSecoes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PopSecoes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in Smooth.
function Smooth_Callback(hObject, eventdata, handles)
% hObject    handle to Smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.SmoothV = get(hObject,'Value');
guidata(hObject,handles)
% save([handles.Directory '\GUIDATA.mat'],'handles')
% Hint: get(hObject,'Value') returns toggle state of Smooth


% --- Executes on button press in Kml2D.
function Kml2D_Callback(hObject, eventdata, handles)
% hObject    handle to Kml2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Kml2DV = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of Kml2D


% --- Executes on button press in Run2D.
function Run2D_Callback(hObject, eventdata, handles)
% hObject    handle to Run2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UserSettings = handles.UserSettings;
UserSettings.Options2D.Smoo2dValue = get(handles.Options2D.Smoo2dOption,'UserData');
UserSettings.Options2D.quiver2dValue = get(handles.Options2D.quiv2dOption,'UserData');
UserSettings.Options2D.plotpostionValue = get(handles.Options2D.plotpostion,'UserData');
UserSettings.Options2D.ExpShp2dValue = get(handles.Options2D.ExpShp2d,'UserData');
UserSettings.Options2D.shpfileexpValue = get(handles.Options2D.shpfileexp,'UserData');
UserSettings.Options2D.extrapnumpointsValue = get(handles.Options2D.extrapnumpoints,'UserData');
UserSettings.Options2D.VortSideValue = get(handles.Options2D.VortSide,'UserData');
UserSettings.Options2D.Z0Ustar1Value = get(handles.Options2D.Z0Ustar,'UserData');
UserSettings.Options2D.NumCellsExtrapValue = get(handles.Options2D.NumCellsExtrap,'UserData');
UserSettings.Options2D.KML2dValue = get(handles.Options2D.KML2doption,'UserData');
p = mfilename('fullpath');
save([p '\UserSettings.mat']);

handles.run2d = 1;
kml2 = handles.Kml2DV;
smoo = handles.SmoothV;
graf = handles.GraficosV;
Sec = handles.PopSecoesV;
exc2 = handles.Excel2DV;
Plot = handles.PlotV;
quiv = handles.quiver2dV;
Directory = handles.Directory;
shp = handles.ShapeV;
zone=[];
if exist([Directory 'Session ' Sec(1:end-4)],'dir')
else
    mkdir([Directory 'Session ' Sec(1:end-4)])
end

if quiv==1
    [Vel1, Vel2]=quiver2DGUI(handles);
    quiv = [{quiv} {Vel1} {Vel2}];
else
    quiv = {quiv};
end
if handles.exp2grid2D.V == 1
    
    [s,v] = listdlg('ListString',[{'Latitude and Longitude'};{'UTM'}],...
        'SelectionMode','single','Name','Grid Input Type','ListSize',[160 35]);
    if v == 0
        return
    end
    if s == 2
        zone = inputdlg('Input UTM Zone (i.e 22 M)','Input',1,{'22 M'});
        if sum(size(zone)==0)~=0
            return
        end        
        zone = zone{1};
    end
    [filenamegrid, pathgrid]= uigetfile({'*.mat'},'pick a file with the desired grid');
    if and(sum(filenamegrid==0),sum(pathgrid==0))
        return
    end

    handles.exp2grid2D.path = pathgrid;
    handles.exp2grid2D.filename = filenamegrid;
    handles.exp2grid2D.zone = zone;
end
[Long, Lat,AverageVelocityVector,excel,handles] = loaddata1graph(Directory,Sec,graf,smoo,Plot,quiv,handles);
if kml2 == 1
    display('Writting KML File...')
    GoogleEarthQuiver([Directory 'Session ' Sec(1:end-4) '\' Sec(1:end-4)],Long,Lat,AverageVelocityVector(:,1),...
                AverageVelocityVector(:,2),AverageVelocityVector(:,3),handles);
    copyfile([pwd '\googleearth\data\redcone.dae'],[Directory 'Session ' Sec(1:end-4) '\redcone.dae'])
    display('Pronto!')
end
if exc2 ==1
    display('Writting Report...')
    path=[Directory 'Session ' Sec(1:end-4) '\excelReport'];
    writeExcelReport2d(excel,path)
    display('Done!')
end
if shp == 1
    display('Writting Shape File...')
    if strcmp(get(handles.Options2D.ExpShp2d,'UserData'),'Geovector')
        FileName = [Directory 'Seção ' Sec(1:end-4) '\ShapeFileGeoVector.shp'];
    else
        FileName = [Directory 'Seção ' Sec(1:end-4) '\ShapeFileGeoPoint.shp'];
    end
    ShpFileWrite(Long,Lat,FileName,handles)
    display('Done')
end
if handles.SideVort.V == 1
    sideVort(handles)
end
% save('C:\Users\Roberta\Desktop\Ricardo\Variables Save\handles.mat','handles')
display('All Done!')


% --- Executes on button press in Batimetria.
function Batimetria_Callback(hObject, eventdata, handles)
% hObject    handle to Batimetria (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.BatimetriaV = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of Batimetria


% --- Executes on button press in Grid.
function Grid_Callback(hObject, eventdata, handles)
% hObject    handle to Grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.GridV = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of Grid


% --- Executes on button press in Superior.
function Superior_Callback(hObject, eventdata, handles)
% hObject    handle to Superior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.SuperiorV = get(hObject,'Value');
guidata(hObject,handles)

% Hint: get(hObject,'Value') returns toggle state of Superior


% --- Executes on button press in quiver3d.
function quiver3d_Callback(hObject, eventdata, handles)
% hObject    handle to quiver3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.quiver3dV = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of quiver3d


% --- Executes on button press in simulacao.
function simulacao_Callback(hObject, eventdata, handles)
% hObject    handle to simulacao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.simulacaoV = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of simulacao


% --- Executes on button press in run3d.
function run3d_Callback(hObject, eventdata, handles)
% hObject    handle to run3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UserSettings = handles.UserSettings;
UserSettings.Options3D.Smoo2dValue = get(handles.Options3D.Smoo3dOption,'UserData');
UserSettings.Options3D.quiver2dValue = get(handles.Options3D.quiv3dOption,'UserData');
UserSettings.Options3D.plotpostionValue = get(handles.Options3D.plotpostion3d,'UserData');
UserSettings.Options3D.ExpShp2dValue = get(handles.Options3D.ExpShp3d,'UserData');
UserSettings.Options3D.shpfileexpValue = get(handles.Options3D.shpfileexp3d,'UserData');
UserSettings.Options3D.Z0Ustar1Value = get(handles.Options3D.Z0Ustar,'UserData');
UserSettings.Options3D.NumCellsExtrapValue = get(handles.Options3D.NumCellsExtrap,'UserData');
UserSettings.Options3D.KML3dValue = get(handles.Options3D.KML3doption,'UserData');
UserSettings.Options3D.extrapnumpointsValue = get(handles.Options3D.quiver3doption,'UserData');
UserSettings.Options3D.smallgraphlimts = get(handles.Options3D.smallgraphlimts,'UserData');
UserSettings.Options3D.SimulOptValue = get(handles.Options3D.SimulOpt,'UserData');
UserSettings.Options3D.VisTopQuivValue = get(handles.Options3D.VisTopQuiv,'UserData');
UserSettings.Options3D.GridTopQuivValue = get(handles.Options3D.GridTopQuiv,'UserData');
UserSettings.Options3D.ContourLinesValue = get(handles.Options3D.ContourLines,'UserData');
UserSettings.Options3D.ContourMapValue = get(handles.Options3D.ContourMap,'UserData');
UserSettings.Options3D.quiver3dValue = get(handles.Options3D.quiver3doption,'UserData');
UserSettings.Options3D.NumofPointsBatValue = get(handles.Options3D.NumofPointsBat,'UserData');
p = mfilename('fullpath');
save([p '\UserSettings.mat']);


handles.run3d = 1;
tables = handles.tabledata;
batimetriav = handles.BatimetriaV;
grid = handles.GridV;
super = handles.SuperiorV;
quiv3 = handles.quiver3dV;
simul = handles.simulacaoV;
kml3 = handles.Kml3dV;
exc3 = handles.Excel3DV;
PlSav = handles.PlotSaveAllV;
Shp3 = handles.Shape3DV;
quiv2D = handles.quiver2d3dV;
smoo = handles.Smooth2d3dV;
StreamLines = handles.StreamLinesV;
plotpreview = handles.plotpreview;
a = tables(:,2);
b= true(size(a));
zone=[];
for i=1:length(a)
    if or(a{i}==1,strcmp(a{i},'true'))
        b(i) = true;
    elseif or(a{i}==0,strcmp(a{i},'false'))
        b(i) = false;
    end
end
if quiv2D==1
    [Vel1, Vel2]=quiver2DGUI(handles);
    quiv2D = [{quiv2D} {Vel1} {Vel2}];
else
    quiv2D = {quiv2D};
end

t= tables(:,1);
OpenFiles = t(b);
Directory = handles.Directory;
if PlSav ==1
    figure
end
if batimetriav == 1
    [filenamebat, pathbat]= uigetfile({'*.xyz'},'Pick a file with bathymetry of river');
    if and(sum(filenamebat==0),sum(pathbat==0))
        return
    end
else
    filenamebat=[];pathbat=[];
end
if grid == 1
    [s,v] = listdlg('ListString',[{'Latitude and Longitude'};{'UTM'}],...
        'SelectionMode','single','Name','Grid Input Type','ListSize',[160 35]);
    if v == 0
        return
    end
    if s == 2
        zone = inputdlg('Input UTM Zone (i.e 22 M)','Input',1,{'22 M'});
        if sum(size(zone)==0)~=0
            return
        end        
        zone = zone{1};
    end
    [filenamegrid, pathgrid]= uigetfile({'*.mat'},'Pick a file with desired grid');
    if and(sum(filenamegrid==0),sum(pathgrid==0))
        return
    end
else
    filenamegrid=[];pathgrid=[];
end
[Quiver3VectorXTot,Quiver3VectorYTot,Quiver3VectorZTot,Quiver3VectorUTot,...
    Quiver3VectorVTot,Quiver3VectorWTot,YrR,XrR,FrR,alt,excel,GridVar,VarNme] =...
    LoadData3D(Directory,batimetriav,OpenFiles,grid,PlSav,super,quiv3,kml3,smoo,simul,...
    quiv2D,filenamegrid, pathgrid,filenamebat, pathbat,zone,StreamLines,handles);
if handles.plotpreview ==1
    for n = 1:length(OpenFiles)
        save([Directory 'Seção ' OpenFiles{n} '/LongitudeLatitudeMagSeçao'],VarNme{n})
    end
end
if simul == 1
    direcao=SimuVelocidades;
    display('Making animation...')
    simulationRiver3D(Quiver3VectorXTot,Quiver3VectorYTot,...
        Quiver3VectorZTot,Quiver3VectorUTot,Quiver3VectorVTot,Quiver3VectorWTot,...
        YrR,XrR,FrR,alt,batimetriav,direcao,handles);
    display('Done!')
end
if exc3 == 1
    display('Writing Report...')
    path=[Directory 'Analysis 3D\ExcelReportAll'];
    writeExcelReport(excel,path)
    display('Done!')
    
end
if Shp3==1
    display('Writting Shape File...')
    if strcmp(get(handles.Options3D.ExpShp3d,'UserData'),'Geovector')        
        FileName = [Directory 'Analysis 3D\All Sessions ShapeFile GeoVector'];
    else
        FileName = [Directory 'Analysis 3D\All Sessions ShapeFile GeoPoints'];
    end
    ShpFileWrite(Quiver3VectorXTot,Quiver3VectorYTot,FileName,handles)
    display('Shape File done!')
end
if StreamLines ==1
    QuiverWithStreamLines(GridVar)
end
display('All done!')
% save('C:\Users\Roberta\Desktop\Ricardo\Variables Save\handles.mat','handles')



% --- Executes on button press in Kml3d.
function Kml3d_Callback(hObject, eventdata, handles)
% hObject    handle to Kml3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Kml3dV = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of Kml3d


% --- Executes during object creation, after setting all properties.
function Secoes3D_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Secoes3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --- Executes when entered data in editable cell(s) in Secoes3D.
function Secoes3D_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to Secoes3D (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
% save('C:\Users\Roberta\Desktop\Ricardo\Variables Save\data.mat','hObject','eventdata','handles')
editData = eventdata.EditData;
indData = eventdata.Indices;
if editData == 0
    editData = 'false';
elseif editData==1
    editData = 'true';
end
tabledata = handles.tabledata;
tabledata{indData(1),2} = editData;
handles.tabledata = tabledata;
set(handles.Secoes3D,'data',handles.tabledata)
guidata(hObject,handles)


% --- Executes on button press in Excel2D.
function Excel2D_Callback(hObject, eventdata, handles)
% hObject    handle to Excel2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Excel2DV = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of Excel2D


% --- Executes on button press in Excel3D.
function Excel3D_Callback(hObject, eventdata, handles)
% hObject    handle to Excel3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Excel3DV = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of Excel3D


% --- Executes on button press in MudarPasta.
function MudarPasta_Callback(hObject, eventdata, handles)
% hObject    handle to MudarPasta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Directory = uigetdir(pwd);
handles.Directory = [Directory '\'];
OpenFiles = DirectoryOpenFiles(Directory);
handles.OpenFiles = OpenFiles;
set(handles.PopSecoes,'String',OpenFiles)
GraficosPPlotar = [{'North'} {'East'} {'Upstream'} {'Magnitude'} {'Error'} {'Transverse'} {'Longitudinal'}];
handles.GraficosPPlotar = GraficosPPlotar;
set(handles.Graficos,'String',GraficosPPlotar)
ColumnFormat={'char','logical'};
ColumnName = {'Sessions','Process'};
% data = {OpenFiles, true(size(OpenFiles))};
% data(:,2) = {'true'};
data1=[];
for i = 1:length(OpenFiles)
    data1 = [data1;{OpenFiles{i}(1:end-4) true}];
end
set(handles.Secoes3D,'data',data1,'ColumnName',ColumnName,'ColumnFormat',ColumnFormat,...
    'ColumnEditable', [false true],'RowName',[],'ColumnWidth',{110, 75});
handles.tabledata = data1;
if length(OpenFiles)<1
    handles.PopSecoesV = '';
else
    handles.PopSecoesV = OpenFiles{1};
end
DisplayDir_Callback(handles.DisplayDir, eventdata, handles)
guidata(hObject,handles)



% --- Executes on button press in PlotSaveAll.
function PlotSaveAll_Callback(hObject, eventdata, handles)
% hObject    handle to PlotSaveAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.PlotSaveAllV = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of PlotSaveAll


% --- Executes on button press in ShapeFile.
function ShapeFile_Callback(hObject, eventdata, handles)
% hObject    handle to ShapeFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ShapeV = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of ShapeFile


% --- Executes on button press in Shp3D.
function Shp3D_Callback(hObject, eventdata, handles)
% hObject    handle to Shp3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Shape3DV = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of Shp3D



function DisplayDir_Callback(hObject, eventdata, handles)
% hObject    handle to DisplayDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Dir = handles.Directory;
set(hObject,'String',Dir)
% set(handles.DisplayDir,'Eneble','off')
% Hints: get(hObject,'String') returns contents of DisplayDir as text
%        str2double(get(hObject,'String')) returns contents of DisplayDir as a double


% --- Executes during object creation, after setting all properties.
function DisplayDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DisplayDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Plot.
function Plot_Callback(hObject, eventdata, handles)
% hObject    handle to Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.PlotV = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of Plot


% --- Executes on button press in quiver2d.
function quiver2d_Callback(hObject, eventdata, handles)
% hObject    handle to quiver2d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.quiver2dV = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of quiver2d


% --- Executes on button press in Smooth3d2d.
function Smooth3d2d_Callback(hObject, eventdata, handles)
% hObject    handle to Smooth3d2d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Smooth2d3dV = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of Smooth3d2d


% --- Executes on button press in quiver2d3d.
function quiver2d3d_Callback(hObject, eventdata, handles)
% hObject    handle to quiver2d3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.quiver2d3dV = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of quiver2d3d


% --- Executes on button press in StreamLines.
function StreamLines_Callback(hObject, eventdata, handles)
% hObject    handle to StreamLines (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.StreamLinesV = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of StreamLines

% function inSmo2d = InputSmooth2d(varargin)
% data = get(varargin{1},'UserData');
% datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
% inSmo2d = inputdlg({'Numero para multiplicar o numero de celulas em x',...
%     'Numero para multiplicar o numero de celulas em y'},'Input',1,{datastr{1:2}});
% if sum(size(inSmo2d)==0)~=0
%     return
% end
% inSmo2d1 = str2num(inSmo2d{1});
% inSmo2d2 = str2num(inSmo2d{2});
% inSmo2d = [inSmo2d1, inSmo2d2];
% set(varargin{1},'UserData',inSmo2d);
% 
% function quiv = Inputquiver2d(varargin)
% data = get(varargin{1},'UserData');
% datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
% quiv = inputdlg('Quiver escala','Input',1,datastr(1));
% if sum(size(quiv)==0)~=0
%     return
% end
% quiv = str2num(quiv{1});
% set(varargin{1},'UserData',quiv);
% 
% function shp2d = InputExpShp2d(varargin)
% source = varargin{1};
% shp2d = get(source,'Label');
% parent = get(source,'Parent');
% switch shp2d
%     case 'Geopoint'
%         set(parent,'UserData','Geopoint');
%     case 'Geovector'
%         set(parent,'UserData','Geovector');        
% end
% 
% function kml = InputKML2d(varargin)
% data = get(varargin{1},'UserData');
% datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
% kml = inputdlg('KML escala Flecha','Input',1,datastr(1));
% if sum(size(kml)==0)~=0
%     return
% end
% kml = str2num(kml{1});
% set(varargin{1},'UserData',kml);
% 
% function plotpos = PlotPos(varargin)
% data = get(varargin{1},'UserData');
% datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
% plp = inputdlg({'Distancia do lado esquerdo do monitor','Distancia de baixo do monitor',...
%     'Largura','Altura'},'Input',1,{datastr{1:4}});
% if sum(size(plp)==0)~=0
%     return
% end
% plp1 = str2num(plp{1});
% plp2 = str2num(plp{2});
% plp3 = str2num(plp{3});
% plp4 = str2num(plp{4});
% plotpos = [plp1 plp2 plp3 plp4];
% set(varargin{1},'UserData',plotpos);
% 
% function expn = extrpnum(varargin)
% data = get(varargin{1},'UserData');
% datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
% expn = inputdlg('Numero de pontos usados em na extrapolaçao para area lateral','Input',1,datastr(1));
% if sum(size(expn)==0)~=0
%     return
% end
% expn = str2num(expn{1});
% set(varargin{1},'UserData',expn);
% 
% function shpfileexp = shpfileexpfun(varargin)
% source = varargin{1};
% shpfileexp = get(source,'Label');
% parent = get(source,'Parent');
% switch shpfileexp
%     case 'utm'
%         set(parent,'UserData','utm');
%     case 'Lat & Long'
%         set(parent,'UserData','Lat & Long');        
% end
% 
% function Chen = cheneq(varargin)
% source = varargin{1};
% Chen = get(source,'Label');
% parent = get(source,'Parent');
% data = get(parent,'UserData');
% datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
% switch Chen
%     case 'Z_0'
%         input = inputdlg('Valor para eq. de Chen','Input',1,datastr(2));
%         if sum(size(input)==0)~=0
%             return
%         end
%         
%         data(2) = str2num(input{1});
%         set(parent,'UserData',data);
%     case 'U*'
%         input = inputdlg('Valor para eq. de Chen','Input',1,datastr(1));
%         if sum(size(input)==0)~=0
%             return
%         end
%         data(1) = str2num(input{1});
%         set(parent,'UserData',data);   
%     case 'n'
%         input = inputdlg('Valor para eq. de Chen','Input',1,datastr(3));
%         if sum(size(input)==0)~=0
%             return
%         end
%         data(3) = str2num(input{1});
%         set(parent,'UserData',data);    
% end
% 
% function numbat = Inputnumbat(varargin)
% data = get(varargin{1},'UserData');
% datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
% numbat = inputdlg('Numero de pontos','Input',1,datastr(1));
% if sum(size(numbat)==0)~=0
%     return
% end
% numbat = str2num(numbat{1});
% set(varargin{1},'UserData',numbat);
% 
% function minipos = minigap(varargin)
% data = get(varargin{1},'UserData');
% datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
% plp = inputdlg({'X lim menor','X lim maior',...
%     'Y lim menor','Y lim maior'},'Input',1,{datastr{1:4}});
% if sum(size(plp)==0)~=0
%     return
% end
% plp1 = str2num(plp{1});
% plp2 = str2num(plp{2});
% plp3 = str2num(plp{3});
% plp4 = str2num(plp{4});
% minipos = [plp1 plp2 plp3 plp4];
% set(varargin{1},'UserData',minipos);


% --- Executes on button press in plotpreview.
function plotpreview_Callback(hObject, eventdata, handles)
% hObject    handle to plotpreview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.plotpreview = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of plotpreview


% --- Executes on button press in plotpreview2d.
function plotpreview2d_Callback(hObject, eventdata, handles)
% hObject    handle to plotpreview2d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.plotpreview2d = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of plotpreview2d


% --- Executes on button press in topVorticity.
function topVorticity_Callback(hObject, eventdata, handles)
% hObject    handle to topVorticity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.topVorticity = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of topVorticity


% --- Executes on button press in SideVort.
function SideVort_Callback(hObject, eventdata, handles)
% hObject    handle to SideVort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.SideVort.V = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of SideVort

% function sideVErtout = SideVert(varargin)
% source = varargin{1};
% sideVErtout = get(source,'Label');
% parent = get(source,'Parent');
% data = get(parent,'UserData');
% datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
% switch sideVErtout
%     case 'Number of Arrows'
%         input = inputdlg('Number of Arrows:','Input',1,datastr(1));
%         if sum(size(input)==0)~=0
%             return
%         end
%         input = str2num(input{1});
%         set(parent,'UserData',[input data(2)]);
%     case 'Percentage of the length'
%         input = inputdlg('Percentage of the Length:','Input',1,datastr(2));
%         if sum(size(input)==0)~=0
%             return
%         end
%         input = str2num(input{1});
%         if or(input<0,input>1)
%             error('Percentage must be bigger than 0 and smaller than 1')
%         end
%         set(parent,'UserData',[data(1) input]);        
% end
% 
% function SimulOptCall(varargin)
% source = varargin{1};
% label = get(source,'Label');
% parent = get(source,'Parent');
% data = get(parent,'UserData');
% datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
% switch label
%     case 'Pause time'
%         input = inputdlg('Pause time:','Input Pause Time',1,datastr(1));
%         if sum(size(input)==0)~=0
%             return
%         end
%         input = str2num(input{1});
%         set(parent,'UserData',[input data(2:end)]);
%     case 'Vector of Rotation'
%         input = inputdlg({'X', 'Y' , 'Z','Deg'},'Rotate',...
%             1,{datastr{2:5}});
%         if sum(size(input)==0)~=0
%             return
%         end
%         input1 = zeros(1,4);
%         for jj = 1:4
%             input1(jj) = str2num(input{jj});
%         end
%         set(parent,'UserData',[data(1) input1]);        
% end

% --- Executes during object creation, after setting all properties.
function SideVort_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SideVort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% c=uicontextmenu;
% set(hObject,'UIContextMenu',c);
% m1 = uimenu(c,'Label','Number of Arrows','Callback',@SideVert);
% m2 = uimenu(c,'Label','Percentage of the length','Callback',@SideVert);


% --- Executes on button press in exp2grid2D.
function exp2grid2D_Callback(hObject, eventdata, handles)
% hObject    handle to exp2grid2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.exp2grid2D.V = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of exp2grid2D


% --- Executes on button press in multisecgui.
function multisecgui_Callback(hObject, eventdata, handles)
% hObject    handle to multisecgui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
multisessionGUI(handles)

% function quiv = contourlines(varargin)
% data = get(varargin{1},'UserData');
% datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
% quiv = inputdlg('Numero de linhas contour','Input',1,datastr(1));
% if sum(size(quiv)==0)~=0
%     return
% end
% quiv = str2num(quiv{1});
% if rem(quiv,1)~=0 || quiv<=0
%     error('Valor deve ser inteiro, positivo e maior que zero')
% end
% set(varargin{1},'UserData',quiv);
% 
% function quiv = contourmap(varargin)
% data = get(varargin{1},'UserData');
% str = {'gray','parula','jet','hsv','hot','cool','spring','summer','autumn','winter','bone','copper',...
%     'pink','lines','colorcube','prism','flag','white'};
% startIndex = regexp(str,data);
% initialv = find(~cellfun(@isempty,startIndex));
% [s,v] = listdlg('PromptString','Selecione um ColorMap:',...
%                 'SelectionMode','single',...
%                 'ListString',str,'InitialValue',initialv);
% if v==0
%     return
% end
% quiv=str{s};
% set(varargin{1},'UserData',quiv);
% 
% 
% function Reset2D(varargin)
% handles = varargin{3};
% p = mfilename('fullpath');
% UserSettings = load([p '\DefautSettings.mat']);
% UserSettings = UserSettings.UserSettings;
% handles.UserSettings.Options2D = UserSettings.Options2D;
% set(handles.Options2D.Smoo2dOption,'UserData',UserSettings.Options2D.Smoo2dValue);
% set(handles.Options2D.quiv2dOption,'UserData',UserSettings.Options2D.quiver2dValue);
% set(handles.Options2D.plotpostion,'UserData',UserSettings.Options2D.plotpostionValue);
% set(handles.Options2D.ExpShp2d,'UserData',UserSettings.Options2D.ExpShp2dValue);
% set(handles.Options2D.shpfileexp,'UserData',UserSettings.Options2D.shpfileexpValue);
% set(handles.Options2D.extrapnumpoints,'UserData',UserSettings.Options2D.extrapnumpointsValue);
% set(handles.Options2D.VortSide,'UserData',UserSettings.Options2D.VortSideValue);
% set(handles.Options2D.Z0Ustar,'UserData',UserSettings.Options2D.Z0Ustar1Value);
% set(handles.Options2D.NumCellsExtrap,'UserData',UserSettings.Options2D.NumCellsExtrapValue);
% set(handles.Options2D.KML2doption,'UserData',UserSettings.Options2D.KML2dValue);
% save([p '\UserSettings.mat'],'UserSettings')
% guidata(varargin{1},handles)
% 
% function Reset3D(hObject, eventdata, handles)
% p = mfilename('fullpath');
% UserSettings = load([p '\DefautSettings.mat']);
% UserSettings = UserSettings.UserSettings;
% handles.UserSettings.Options3D = UserSettings.Options3D;
% 
% set(handles.Options3D.Smoo3dOption,'UserData',UserSettings.Options3D.Smoo2dValue);
% set(handles.Options3D.quiv3dOption,'UserData',UserSettings.Options3D.quiver2dValue);
% set(handles.Options3D.plotpostion3d,'UserData',UserSettings.Options3D.plotpostionValue);
% set(handles.Options3D.ExpShp3d,'UserData',UserSettings.Options3D.ExpShp2dValue);
% set(handles.Options3D.shpfileexp3d,'UserData',UserSettings.Options3D.shpfileexpValue);
% set(handles.Options3D.Z0Ustar,'UserData',UserSettings.Options3D.Z0Ustar1Value);
% set(handles.Options3D.NumCellsExtrap,'UserData',UserSettings.Options3D.NumCellsExtrapValue);
% set(handles.Options3D.KML3doption,'UserData',UserSettings.Options3D.KML3dValue);
% set(handles.Options3D.quiver3doption,'UserData',UserSettings.Options3D.extrapnumpointsValue);
% set(handles.Options3D.smallgraphlimts,'UserData',UserSettings.Options3D.smallgraphlimts);
% set(handles.Options3D.SimulOpt,'UserData',UserSettings.Options3D.SimulOptValue);
% set(handles.Options3D.VisTopQuiv,'UserData',UserSettings.Options3D.VisTopQuivValue);
% set(handles.Options3D.GridTopQuiv,'UserData',UserSettings.Options3D.GridTopQuivValue);
% set(handles.Options3D.ContourLines,'UserData',UserSettings.Options3D.ContourLinesValue);
% set(handles.Options3D.ContourMap,'UserData',UserSettings.Options3D.ContourMapValue);
% set(handles.Options3D.quiver3doption,'UserData',UserSettings.Options3D.quiver3dValue);
% set(handles.Options3D.NumofPointsBat,'UserData',UserSettings.Options3D.NumofPointsBatValue);
% save([p '\UserSettings.mat'],'UserSettings')
% guidata(hObject,handles)


% --------------------------------------------------------------------
function menu_Callback(hObject, eventdata, handles)
% hObject    handle to menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu1_Callback(hObject, eventdata, handles)
% hObject    handle to menu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function SmoothMenu_Callback(hObject, eventdata, handles)
% hObject    handle to SmoothMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA),
data = get(handles.Options2D.Smoo2dOption,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
inSmo2d = inputdlg({'Number to multiply the number of cells in x',...
    'Number to multiply the number of cells in y'},'Input',1,{datastr{1:2}});
if sum(size(inSmo2d)==0)~=0
    return
end

inSmo2d1 = str2num(inSmo2d{1});
inSmo2d2 = str2num(inSmo2d{2});
inSmo2d = [inSmo2d1, inSmo2d2];
if rem(inSmo2d1,1)~=0 || inSmo2d1==0
    error('Must be an integer greater than zero (x>0)')
end
if rem(inSmo2d2,1)~=0 || inSmo2d2==0
    error('Must be an integer greater than zero (x>0)')
end
set(handles.Options2D.Smoo2dOption,'UserData',inSmo2d);



% --------------------------------------------------------------------
function QuiverMenu_Callback(hObject, eventdata, handles)
% hObject    handle to QuiverMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.Options2D.quiv2dOption,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
quiv = inputdlg('Quiver scale','Input',1,datastr(1));
if sum(size(quiv)==0)~=0
    return
end
quiv = str2num(quiv{1});
set(handles.Options2D.quiv2dOption,'UserData',quiv);

% --------------------------------------------------------------------
function PlotPositionMenu_Callback(hObject, eventdata, handles)
% hObject    handle to PlotPositionMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.Options2D.plotpostion,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
plp = inputdlg({'Left distance from beggining of monitor','Distance from bottom of monitor',...
    'Length','Height'},'Input',1,{datastr{1:4}});
if sum(size(plp)==0)~=0
    return
end
plp1=str2num(plp{1});
plp2=str2num(plp{2});
plp3=str2num(plp{3});
plp4=str2num(plp{4});
% if plp1>plp2
%     xlim = [plp2 plp1];
% else
%     xlim = [plp1 plp2];
% end
% if plp3>plp4
%     ylim = [plp4 plp3];
% else
%     ylim = [plp3 plp4];
% end
plotpos = [plp1 plp2 plp3 plp4];
set(handles.Options2D.plotpostion,'UserData',plotpos);
guidata(hObject,handles)

% --------------------------------------------------------------------
function SideExtrapMenu_Callback(hObject, eventdata, handles)
% hObject    handle to SideExtrapMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.Options2D.extrapnumpoints,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
expn = inputdlg('Number of points for side extrapolation','Input',1,datastr(1));
if sum(size(expn)==0)~=0
    return
end
expn = str2num(expn{1});
if rem(expn,1)~=0 || expn==0
    error('Value must be integer greater than zero')    
end
set(handles.Options2D.extrapnumpoints,'UserData',expn);

% --------------------------------------------------------------------
function VelExtrapMenu_Callback(hObject, eventdata, handles)
% hObject    handle to VelExtrapMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.Options2D.NumCellsExtrap,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
numbat = inputdlg('Number of Points','Input',1,datastr(1));
if sum(size(numbat)==0)~=0
    return
end
numbat = str2num(numbat{1});
if rem(numbat,1)~=0 || numbat==0
    error('Value must be integer greater than zero')    
end
set(handles.Options2D.NumCellsExtrap,'UserData',numbat);
guidata(hObject,handles)

% --------------------------------------------------------------------
function ChenMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ChenMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function SideViewMenu_Callback(hObject, eventdata, handles)
% hObject    handle to SideViewMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function KmlMenu_Callback(hObject, eventdata, handles)
% hObject    handle to KmlMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.Options2D.KML2doption,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
kml = inputdlg('KML arrow scale','Input',1,datastr(1));
if sum(size(kml)==0)~=0
    return
end
kml = str2num(kml{1});
set(handles.Options2D.KML2doption,'UserData',kml);
guidata(hObject,handles)

% --------------------------------------------------------------------
function ShapeMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ShapeMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function ResetDefault2D_Callback(hObject, eventdata, handles)
% hObject    handle to ResetDefault2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
p = mfilename('fullpath');
UserSettings = load([p '\DefautSettings.mat']);
UserSettings = UserSettings.UserSettings;
handles.UserSettings.Options2D = UserSettings.Options2D;
set(handles.Options2D.Smoo2dOption,'UserData',UserSettings.Options2D.Smoo2dValue);
set(handles.Options2D.quiv2dOption,'UserData',UserSettings.Options2D.quiver2dValue);
set(handles.Options2D.plotpostion,'UserData',UserSettings.Options2D.plotpostionValue);
set(handles.Options2D.ExpShp2d,'UserData',UserSettings.Options2D.ExpShp2dValue);
set(handles.Options2D.shpfileexp,'UserData',UserSettings.Options2D.shpfileexpValue);
set(handles.Options2D.extrapnumpoints,'UserData',UserSettings.Options2D.extrapnumpointsValue);
set(handles.Options2D.VortSide,'UserData',UserSettings.Options2D.VortSideValue);
set(handles.Options2D.Z0Ustar,'UserData',UserSettings.Options2D.Z0Ustar1Value);
set(handles.Options2D.NumCellsExtrap,'UserData',UserSettings.Options2D.NumCellsExtrapValue);
set(handles.Options2D.KML2doption,'UserData',UserSettings.Options2D.KML2dValue);
save([p '\UserSettings.mat'],'UserSettings')
guidata(hObject,handles)


% --------------------------------------------------------------------
function ShapeType_Callback(hObject, eventdata, handles)
% hObject    handle to ShapeType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ShpFileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ShpFileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function UstarMenu_Callback(hObject, eventdata, handles)
% hObject    handle to UstarMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.Options2D.Z0Ustar,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
input = inputdlg('Value for U* in Chens equation','Input',1,datastr(1));
if sum(size(input)==0)~=0
    return
end
data(1) = str2num(input{1});
set(handles.Options2D.Z0Ustar,'UserData',data);
guidata(hObject,handles)


% --------------------------------------------------------------------
function ZzeroMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ZzeroMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.Options2D.Z0Ustar,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
input = inputdlg('Value for Z_0 in Chens equation','Input',1,datastr(2));
if sum(size(input)==0)~=0
    return
end
data(2) = str2num(input{1});
set(handles.Options2D.Z0Ustar,'UserData',data);
guidata(hObject,handles)


% --------------------------------------------------------------------
function nChenMenu_Callback(hObject, eventdata, handles)
% hObject    handle to nChenMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.Options2D.Z0Ustar,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
input = inputdlg('Value for n in Chens equation','Input',1,datastr(3));
if sum(size(input)==0)~=0
    return
end
data(3) = str2num(input{1});
set(handles.Options2D.Z0Ustar,'UserData',data);
guidata(hObject,handles)


% --------------------------------------------------------------------
function utmShpOut_Callback(hObject, eventdata, handles)
% hObject    handle to utmShpOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.Options2D.shpfileexp,'UserData','utm')
guidata(hObject,handles)

% --------------------------------------------------------------------
function LongLatShpOut_Callback(hObject, eventdata, handles)
% hObject    handle to LongLatShpOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.Options2D.shpfileexp,'UserData','Lat & Long')
guidata(hObject,handles)


% --------------------------------------------------------------------
function GeopointMenu_Callback(hObject, eventdata, handles)
% hObject    handle to GeopointMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.Options2D.ExpShp2d,'UserData','Geopoint')
guidata(hObject,handles)


% --------------------------------------------------------------------
function GeovectorMenu_Callback(hObject, eventdata, handles)
% hObject    handle to GeovectorMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.Options2D.ExpShp2d,'UserData','Geovector')
guidata(hObject,handles)



% --------------------------------------------------------------------
function Smooth2DMenu3D_Callback(hObject, eventdata, handles)
% hObject    handle to Smooth2DMenu3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.Options3D.Smoo3dOption,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
inSmo2d = inputdlg({'Number to multiply the number of cells in x',...
    'Number to multiply the number of cells in y'},'Input',1,{datastr{1:2}});
if sum(size(inSmo2d)==0)~=0
    return
end

inSmo2d1 = str2num(inSmo2d{1});
inSmo2d2 = str2num(inSmo2d{2});
inSmo2d = [inSmo2d1, inSmo2d2];
if rem(inSmo2d1,1)~=0 || inSmo2d1==0
    error('Must be an integer greater than zero (x>0)')
end
if rem(inSmo2d2,1)~=0 || inSmo2d2==0
    error('Must be an integer greater than zero (x>0)')
end
set(handles.Options3D.Smoo3dOption,'UserData',inSmo2d);


% --------------------------------------------------------------------
function Quiver2DMenu3D_Callback(hObject, eventdata, handles)
% hObject    handle to Quiver2DMenu3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.Options3D.quiv3dOption,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
quiv = inputdlg('Quiver scale','Input',1,datastr(1));
if sum(size(quiv)==0)~=0
    return
end
quiv = str2num(quiv{1});
set(handles.Options3D.quiv3dOption,'UserData',quiv);
guidata(hObject,handles)


% --------------------------------------------------------------------
function ShpMenu3D_Callback(hObject, eventdata, handles)
% hObject    handle to ShpMenu3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function KmlMenu3D_Callback(hObject, eventdata, handles)
% hObject    handle to KmlMenu3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.Options3D.KML3doption,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
kml = inputdlg('KML arrow scale','Input',1,datastr(1));
if sum(size(kml)==0)~=0
    return
end
kml = str2num(kml{1});
set(handles.Options3D.KML3doption,'UserData',kml);
guidata(hObject,handles)


% --------------------------------------------------------------------
function Quiver3DMenu3D_Callback(hObject, eventdata, handles)
% hObject    handle to Quiver3DMenu3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.Options3D.quiver3doption,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
quiv = inputdlg('Quiver scale','Input',1,datastr(1));
if sum(size(quiv)==0)~=0
    return
end
quiv = str2num(quiv{1});
set(handles.Options3D.quiver3doption,'UserData',quiv);
guidata(hObject,handles)

% --------------------------------------------------------------------
function BathPointsMenu3D_Callback(hObject, eventdata, handles)
% hObject    handle to BathPointsMenu3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.Options3D.NumofPointsBat,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
numbat = inputdlg('Number of points','Input',1,datastr(1));
if sum(size(numbat)==0)~=0
    return
end
numbat = str2num(numbat{1});
if rem(numbat,1)~=0 || numbat==0
    error('Must be an integer greater than zero (x>0)')
end
set(handles.Options3D.NumofPointsBat,'UserData',numbat);
guidata(hObject,handles)



% --------------------------------------------------------------------
function MiniPLotMenu3D_Callback(hObject, eventdata, handles)
% hObject    handle to MiniPLotMenu3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.Options3D.smallgraphlimts,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
plp = inputdlg({'X lim small','X lim big',...
    'Y lim small','Y lim big'},'Input',1,{datastr{1:4}});
if sum(size(plp)==0)~=0
    return
end
plp1 = str2num(plp{1});
plp2 = str2num(plp{2});
plp3 = str2num(plp{3});
plp4 = str2num(plp{4});
if plp1>plp2
    xlim = [plp2 plp1];
else
    xlim = [plp1 plp2];
end
if plp3>plp4
    ylim = [plp4 plp3];
else
    ylim = [plp3 plp4];
end
minipos = [xlim ylim];
set(handles.Options3D.smallgraphlimts,'UserData',minipos);
guidata(hObject,handles)

% --------------------------------------------------------------------
function PlotPositionMenu3D_Callback(hObject, eventdata, handles)
% hObject    handle to PlotPositionMenu3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.Options3D.plotpostion3d,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
plp = inputdlg({'Left distance from beggining of monitor','Distance from bottom of monitor',...
    'Length','Height'},'Input',1,{datastr{1:4}});
if sum(size(plp)==0)~=0
    return
end
plp1=str2num(plp{1});
plp2=str2num(plp{2});
plp3=str2num(plp{3});
plp4=str2num(plp{4});
% if plp1>plp2
%     xlim = [plp2 plp1];
% else
%     xlim = [plp1 plp2];
% end
% if plp3>plp4
%     ylim = [plp4 plp3];
% else
%     ylim = [plp3 plp4];
% end
plotpos = [plp1 plp2 plp3 plp4];
set(handles.Options3D.plotpostion3d,'UserData',plotpos);
guidata(hObject,handles)


% --------------------------------------------------------------------
function VelExtrapMenu3D_Callback(hObject, eventdata, handles)
% hObject    handle to VelExtrapMenu3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.Options3D.NumCellsExtrap,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
numbat = inputdlg('Number of Points','Input',1,datastr(1));
if sum(size(numbat)==0)~=0
    return
end
numbat = str2num(numbat{1});
if rem(numbat,1)~=0 || numbat==0
    error('Value must be integer greater than zero')    
end
set(handles.Options3D.NumCellsExtrap,'UserData',numbat);
guidata(hObject,handles)


% --------------------------------------------------------------------
function ChenMenu3D_Callback(hObject, eventdata, handles)
% hObject    handle to ChenMenu3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function AnimationMenu3D_Callback(hObject, eventdata, handles)
% hObject    handle to AnimationMenu3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function TopViewQuivMenu3D_Callback(hObject, eventdata, handles)
% hObject    handle to TopViewQuivMenu3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.Options3D.VisTopQuiv,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
quiv = inputdlg('Quiver scale','Input',1,datastr(1));
if sum(size(quiv)==0)~=0
    return
end
quiv = str2num(quiv{1});
set(handles.Options3D.VisTopQuiv,'UserData',quiv);
guidata(hObject,handles)


% --------------------------------------------------------------------
function TopViewGridMenu3D_Callback(hObject, eventdata, handles)
% hObject    handle to TopViewGridMenu3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.Options3D.GridTopQuiv,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
quiv = inputdlg('Quiver scale','Input',1,datastr(1));
if sum(size(quiv)==0)~=0
    return
end
quiv = str2num(quiv{1});
set(handles.Options3D.GridTopQuiv,'UserData',quiv);
guidata(hObject,handles)

% --------------------------------------------------------------------
function NumColorBathMenu3D_Callback(hObject, eventdata, handles)
% hObject    handle to NumColorBathMenu3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.Options3D.ContourLines,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
quiv = inputdlg('Number of lines in contour of bathymetry','Input',1,datastr(1));
if sum(size(quiv)==0)~=0
    return
end
quiv = str2num(quiv{1});
if rem(quiv,1)~=0 || quiv<=0
    error('Must be an integer greater than zero (x>0)')
end
set(handles.Options3D.ContourLines,'UserData',quiv);
guidata(hObject,handles)



% --------------------------------------------------------------------
function ColorMapMenu3D_Callback(hObject, eventdata, handles)
% hObject    handle to ColorMapMenu3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.Options3D.ContourMap,'UserData');
str = {'gray','parula','jet','hsv','hot','cool','spring','summer','autumn','winter','bone','copper',...
    'pink','lines','colorcube','prism','flag','white'};
startIndex = regexp(str,data);
initialv = find(~cellfun(@isempty,startIndex));
[s,v] = listdlg('PromptString','Select a ColorMap:',...
                'SelectionMode','single',...
                'ListString',str,'InitialValue',initialv);
if v==0
    return
end
quiv=str{s};
set(handles.Options3D.ContourMap,'UserData',quiv);

% --------------------------------------------------------------------
function Reset3D_Callback(hObject, eventdata, handles)
% hObject    handle to Reset3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
p = mfilename('fullpath');
UserSettings = load([p '\DefautSettings.mat']);
UserSettings = UserSettings.UserSettings;
handles.UserSettings.Options3D = UserSettings.Options3D;

set(handles.Options3D.Smoo3dOption,'UserData',UserSettings.Options3D.Smoo2dValue);
set(handles.Options3D.quiv3dOption,'UserData',UserSettings.Options3D.quiver2dValue);
set(handles.Options3D.plotpostion3d,'UserData',UserSettings.Options3D.plotpostionValue);
set(handles.Options3D.ExpShp3d,'UserData',UserSettings.Options3D.ExpShp2dValue);
set(handles.Options3D.shpfileexp3d,'UserData',UserSettings.Options3D.shpfileexpValue);
set(handles.Options3D.Z0Ustar,'UserData',UserSettings.Options3D.Z0Ustar1Value);
set(handles.Options3D.NumCellsExtrap,'UserData',UserSettings.Options3D.NumCellsExtrapValue);
set(handles.Options3D.KML3doption,'UserData',UserSettings.Options3D.KML3dValue);
set(handles.Options3D.quiver3doption,'UserData',UserSettings.Options3D.extrapnumpointsValue);
set(handles.Options3D.smallgraphlimts,'UserData',UserSettings.Options3D.smallgraphlimts);
set(handles.Options3D.SimulOpt,'UserData',UserSettings.Options3D.SimulOptValue);
set(handles.Options3D.VisTopQuiv,'UserData',UserSettings.Options3D.VisTopQuivValue);
set(handles.Options3D.GridTopQuiv,'UserData',UserSettings.Options3D.GridTopQuivValue);
set(handles.Options3D.ContourLines,'UserData',UserSettings.Options3D.ContourLinesValue);
set(handles.Options3D.ContourMap,'UserData',UserSettings.Options3D.ContourMapValue);
set(handles.Options3D.quiver3doption,'UserData',UserSettings.Options3D.quiver3dValue);
set(handles.Options3D.NumofPointsBat,'UserData',UserSettings.Options3D.NumofPointsBatValue);
save([p '\UserSettings.mat'],'UserSettings')
guidata(hObject,handles)


% --------------------------------------------------------------------
function PauseAniMenu3D_Callback(hObject, eventdata, handles)
% hObject    handle to PauseAniMenu3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.Options3D.SimulOpt,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
input = inputdlg('Pause time:','Input Pause Time',1,datastr(1));
if sum(size(input)==0)~=0
    return
end
input = str2num(input{1});
set(handles.Options3D.SimulOpt,'UserData',[input data(2:end)]); 
guidata(hObject,handles)



% --------------------------------------------------------------------
function RotAniMenu3D_Callback(hObject, eventdata, handles)
% hObject    handle to RotAniMenu3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = get(handles.Options3D.SimulOpt,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
input = inputdlg({'X', 'Y' , 'Z','Deg'},'Rotate',...
    1,datastr(2:5));
if sum(size(input)==0)~=0
    return
end
input1 = zeros(1,4);
for jj = 1:4
    input1(jj) = str2num(input{jj});
end
set(handles.Options3D.SimulOpt,'UserData',[data(1) input1]); 
guidata(hObject,handles)



% --------------------------------------------------------------------
function UstarMenu3D_Callback(hObject, eventdata, handles)
% hObject    handle to UstarMenu3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.Options3D.Z0Ustar,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
input = inputdlg('Value for U* in Chens equation','Input',1,datastr(1));
if sum(size(input)==0)~=0
    return
end
data(1) = str2num(input{1});
set(handles.Options3D.Z0Ustar,'UserData',data);
guidata(hObject,handles)


% --------------------------------------------------------------------
function ZzeroMenu3D_Callback(hObject, eventdata, handles)
% hObject    handle to ZzeroMenu3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.Options3D.Z0Ustar,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
input = inputdlg('Value for Z_0 in Chens equation','Input',1,datastr(2));
if sum(size(input)==0)~=0
    return
end
data(2) = str2num(input{1});
set(handles.Options3D.Z0Ustar,'UserData',data);
guidata(hObject,handles)


% --------------------------------------------------------------------
function nChenMenu3D_Callback(hObject, eventdata, handles)
% hObject    handle to nChenMenu3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.Options3D.Z0Ustar,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
input = inputdlg('Value for n in Chens equation','Input',1,datastr(3));
if sum(size(input)==0)~=0
    return
end
data(3) = str2num(input{1});
set(handles.Options3D.Z0Ustar,'UserData',data);
guidata(hObject,handles)



% --------------------------------------------------------------------
function TypeShpMenu3D_Callback(hObject, eventdata, handles)
% hObject    handle to TypeShpMenu3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function FileShpMenu3D_Callback(hObject, eventdata, handles)
% hObject    handle to FileShpMenu3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function utmShpMenu3D_Callback(hObject, eventdata, handles)
% hObject    handle to utmShpMenu3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.Options3D.shpfileexp3d,'UserData','utm')
guidata(hObject,handles)

% --------------------------------------------------------------------
function LongLatShpMenu3D_Callback(hObject, eventdata, handles)
% hObject    handle to LongLatShpMenu3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.Options3D.shpfileexp3d,'UserData','Lat & Long')
guidata(hObject,handles)

% --------------------------------------------------------------------
function ShpGeopointMenu3D_Callback(hObject, eventdata, handles)
% hObject    handle to ShpGeopointMenu3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.Options3D.ExpShp3d,'UserData','Geopoint')
guidata(hObject,handles)


% --------------------------------------------------------------------
function ShpGeovectorMenu3D_Callback(hObject, eventdata, handles)
% hObject    handle to ShpGeovectorMenu3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.Options3D.ExpShp3d,'UserData','Geovector')
guidata(hObject,handles)


% --------------------------------------------------------------------
function NumArrowsVortMenu_Callback(hObject, eventdata, handles)
% hObject    handle to NumArrowsVortMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.Options2D.VortSide,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
input = inputdlg('Number of Arrows:','Input',1,datastr(1));
if sum(size(input)==0)~=0
    return
end
input = str2num(input{1});
if rem(input,1)~=0 || input==0
    error('Value must be integer greater than zero')    
end
set(handles.Options2D.VortSide,'UserData',[input data(2)]);


% --------------------------------------------------------------------
function PercentLengthVortMenu_Callback(hObject, eventdata, handles)
% hObject    handle to PercentLengthVortMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.Options2D.VortSide,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
input = inputdlg('Percentage of the Length:','Input',1,datastr(2));
if sum(size(input)==0)~=0
    return
end
input = str2num(input{1});
if or(input<0,input>1)
    error('Percentage must be bigger than 0 and smaller than 1')
end
set(handles.Options2D.VortSide,'UserData',[data(1) input]);


% --------------------------------------------------------------------
function ContactTagUiContext_Callback(hObject, eventdata, handles)
% hObject    handle to ContactTagUiContext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('Contact Ricardo Hopker through: rbhopker@broncs.utpa.edu','Contact')


% --------------------------------------------------------------------
function ContexMenuRicardoHopker_Callback(hObject, eventdata, handles)
% hObject    handle to ContexMenuRicardoHopker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
