function varargout = multisessionGUI(varargin)
% MULTISESSIONGUI MATLAB code for multisessionGUI.fig
%      MULTISESSIONGUI, by itself, creates a new MULTISESSIONGUI or raises the existing
%      singleton*.
%
%      H = MULTISESSIONGUI returns the handle to a new MULTISESSIONGUI or the handle to
%      the existing singleton*.
%
%      MULTISESSIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MULTISESSIONGUI.M with the given input arguments.
%
%      MULTISESSIONGUI('Property','Value',...) creates a new MULTISESSIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before multisessionGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to multisessionGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help multisessionGUI

% Last Modified by GUIDE v2.5 18-Jun-2015 16:04:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @multisessionGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @multisessionGUI_OutputFcn, ...
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


% --- Executes just before multisessionGUI is made visible.
function multisessionGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to multisessionGUI (see VARARGIN)

% Choose default command line output for multisessionGUI
handles.output = hObject;
handles.input = varargin{1};
Directory = handles.input.Directory;
OpenFiles = DirectoryOpenFiles(Directory);
handles.Directory = Directory;
handles.OpenFiles = OpenFiles;
ColumnFormat={'char','logical'};
ColumnName = {'Sections','Process'};
data1=[];
menu = uimenu(gcf,'Label','Options');
NumPointsX = uimenu(menu,'Label','Number of cells through the section','UserData',[100],'Callback',@InputNum);
handles.options.NumPointsX = NumPointsX;
Extrap = uimenu(menu,'Label','Number of cell for side extrapolation','UserData',[10],'Callback',@InputNum);
handles.Options2D.extrapnumpoints = Extrap;
Z0Ustar = uimenu(menu,'Label','Chens equation, Extrapolation','UserData',[.2 .1 1/6]);
Ustar = uimenu(Z0Ustar,'Label','U*','Callback',@cheneq);
Z0 = uimenu(Z0Ustar,'Label','Z_0','Callback',@cheneq);
power = uimenu(Z0Ustar,'Label','n','Callback',@cheneq);
handles.Options3D.Z0Ustar = Z0Ustar;
NumCellsExtrap = uimenu(menu,'Label','Number of cell for velocity extrapolation',...
    'UserData',[5],'Callback',@Inputnumbat);
handles.Options3D.NumCellsExtrap = NumCellsExtrap;
QuivMulti = uimenu(menu,'Label','Number for quiver scale',...
    'UserData',[1],'Callback',@Inputnumbat);
handles.Options3D.QuivMulti = QuivMulti;
KML2doption = uimenu(menu,'Label','Number for KML arrow Scale',...
    'UserData',[100],'Callback',@Inputnumbat);
handles.Options2D.KML2doption = KML2doption;
handles.run3d =0;
shpfileexp = uimenu(menu,'Label','Shape File Export Format','UserData','utm');
shpfileexp1 = uimenu(shpfileexp,'Label','utm','Callback',@shpfileexpfun);
shpfileexp2 = uimenu(shpfileexp,'Label','Lat & Long','Callback',@shpfileexpfun);
ExpShp2d = uimenu(menu,'Label','Shape','UserData','Geovector');
ExpShp2doption1 = uimenu(ExpShp2d,'Label','Geopoint','Callback',@InputExpShp2d);
ExpShp2doption2 = uimenu(ExpShp2d,'Label','Geovector','Callback',@InputExpShp2d);
plotpostion = uimenu(menu,'Label','Plot Position','UserData',[200, 558, 560*3, 420*.75],'Callback',...
    @PlotPos);
handles.Options2D.plotpostion = plotpostion;
VortSide = uimenu(menu,'Label','Side View Vorticity','UserData',[10 0.05]);
VortSide1 = uimenu(VortSide,'Label','Number of Arrows','Callback',@SideVert);
VortSide2 = uimenu(VortSide,'Label','Percentage of the length','Callback',@SideVert);
handles.Options2D.VortSide =VortSide;

handles.Options2D.shpfileexp = shpfileexp;
handles.Options2D.ExpShp2d = ExpShp2d;
for i = 1:length(OpenFiles)
    data1 = [data1;{OpenFiles{i}(1:end-4) true}];
end
set(handles.uitable1,'data',data1,'ColumnName',ColumnName,'ColumnFormat',ColumnFormat,...
    'ColumnEditable', [false true],'RowName',[],'ColumnWidth',{110, 75});
handles.tabledata = data1;
handles.VortLate.V = 0;
handles.GridExp.V = 0;
handles.ShpExp.V = 0;
handles.kmlexp.V = 0;
handles.excelgen.V = 0;
handles.Quiverplot.V = 0;
handles.plotsess.V = 1;

GraficosPPlotar = [{'North'} {'East'} {'Upstream'} {'Magnitude'} {'Error'} {'Transverse'} {'Longitudinal'}];
handles.GraficosPPlotar = GraficosPPlotar;
set(handles.VelPopUp,'String',GraficosPPlotar)
handles.VelPopUp.V = 'North';
displaydir_Callback(handles.displaydir, eventdata, handles)
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes multisessionGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = multisessionGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in ChangeFold.
function ChangeFold_Callback(hObject, eventdata, handles)
% hObject    handle to ChangeFold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Directory = uigetdir(pwd);
handles.Directory = [Directory '\'];
OpenFiles = DirectoryOpenFiles(Directory);
handles.OpenFiles = OpenFiles;
displaydir_Callback(handles.displaydir, eventdata, handles)

ColumnFormat={'char','logical'};
ColumnName = {'Sections','Process'};
data1=[];
for i = 1:length(OpenFiles)
    data1 = [data1;{OpenFiles{i}(1:end-4) true}];
end
set(handles.uitable1,'data',data1,'ColumnName',ColumnName,'ColumnFormat',ColumnFormat,...
    'ColumnEditable', [false true],'RowName',[],'ColumnWidth',{110, 75});
handles.tabledata = data1;

guidata(hObject,handles)


% --- Executes on button press in replot.
function replot_Callback(hObject, eventdata, handles)
% hObject    handle to replot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tables = handles.tabledata;
a = tables(:,2);
b= true(size(a));
for i=1:length(a)
    if or(a{i}==1,strcmp(a{i},'true'))
        b(i) = true;
    elseif or(a{i}==0,strcmp(a{i},'false'))
        b(i) = false;
    end
end
t= tables(:,1);
OpenFiles = t(b);
handles.OpenFiles = OpenFiles;
replot = replot4gui(handles);
Long = replot.Long;
Lat = replot.Lat;
QuiverX = replot.Quiver3X;
QuiverY = replot.Quiver3Y;
QuiverU = replot.Quiver3U;
QuiverV = replot.Quiver3V;
StartEdge = replot.StartEdge;
colours = distinguishable_colors(length(Long)+1);
legendstr = [];
r= zeros(length(Long),1);
mineachX = r; mineachY = r; maxeachX = r; maxeachY = r;
cla(handles.axes1)
for q = 1:length(Long)
    plot(handles.axes1,Long{q},Lat{q},'Color',colours(q,:));
    hold on
    r(q)= quiver(handles.axes1,QuiverX{q},QuiverY{q},QuiverU{q},QuiverV{q},'Color',colours(q,:));
    if StartEdge{q} == 0
        mineachX(q) = QuiverX{q}(1);
        mineachY(q) = QuiverY{q}(1);
        maxeachX(q) = QuiverX{q}(end);
        maxeachY(q) = QuiverY{q}(end);
    else
        mineachX(q) = QuiverX{q}(end);
        mineachY(q) = QuiverY{q}(end);
        maxeachX(q) = QuiverX{q}(1);
        maxeachY(q) = QuiverY{q}(1);
    end
    set(r(q),'ShowArrowHead','off')
    legendstr = [legendstr OpenFiles(q)];  
end
NumOfPoints = get(handles.options.NumPointsX,'UserData');
aveXmin = mean(mineachX);
aveYmin = mean(mineachY);
aveXmax = mean(maxeachX);
aveYmax = mean(maxeachY);
X = linspace(aveXmin,aveXmax,NumOfPoints);
Y = linspace(aveYmin,aveYmax,NumOfPoints);
r(q+1)=plot(handles.axes1,X,Y,'Color',colours(q+1,:));
legendstr = [legendstr {'Mean Path'}];
legend(r,legendstr)
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

% n=get(handles.axes1,'xtick');
% set(handles.axes1,'xticklabel',sprintf('%.4f |',n'));



% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
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
set(handles.uitable1,'data',handles.tabledata)
guidata(hObject,handles)



function displaydir_Callback(hObject, eventdata, handles)
% hObject    handle to displaydir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Dir = handles.Directory;
set(hObject,'String',Dir)

% Hints: get(hObject,'String') returns contents of displaydir as text
%        str2double(get(hObject,'String')) returns contents of displaydir as a double


% --- Executes during object creation, after setting all properties.
function displaydir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to displaydir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function InputNum(varargin)
data = get(varargin{1},'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
numbat = inputdlg('Number of points','Input',1,datastr(1));
if sum(size(numbat)==0)~=0
    return
end
numbat = str2num(numbat{1});
set(varargin{1},'UserData',numbat);


% --- Executes on button press in plotsess.
function plotsess_Callback(hObject, eventdata, handles)
% hObject    handle to plotsess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.plotsess.V = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of plotsess


% --- Executes on button press in Quiverplot.
function Quiverplot_Callback(hObject, eventdata, handles)
% hObject    handle to Quiverplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Quiverplot.V = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of Quiverplot


% --- Executes on button press in excelgen.
function excelgen_Callback(hObject, eventdata, handles)
% hObject    handle to excelgen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.excelgen.V = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of excelgen


% --- Executes on button press in kmlexp.
function kmlexp_Callback(hObject, eventdata, handles)
% hObject    handle to kmlexp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.kmlexp.V = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of kmlexp


% --- Executes on button press in ShpExp.
function ShpExp_Callback(hObject, eventdata, handles)
% hObject    handle to ShpExp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ShpExp.V = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of ShpExp


% --- Executes on button press in GridExp.
function GridExp_Callback(hObject, eventdata, handles)
% hObject    handle to GridExp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.GridExp.V = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of GridExp


% --- Executes on button press in VortLate.
function VortLate_Callback(hObject, eventdata, handles)
% hObject    handle to VortLate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.VortLate.V = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of VortLate


% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tables = handles.tabledata;
a = tables(:,2);
b= true(size(a));
for i=1:length(a)
    if or(a{i}==1,strcmp(a{i},'true'))
        b(i) = true;
    elseif or(a{i}==0,strcmp(a{i},'false'))
        b(i) = false;
    end
end
t= tables(:,1);
OpenFiles = t(b);
handles.OpenFiles = OpenFiles;
quiv2D = handles.Quiverplot.V;
if quiv2D==1
    [Vel1, Vel2]=quiver2DGUI(handles);
    quiv2D = [{quiv2D} {Vel1} {Vel2}];
else
    quiv2D = {quiv2D};
end
handles.Quiverplot.quiv2D = quiv2D;

if handles.GridExp.V == 1
    
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

    handles.GridExp.path = pathgrid;
    handles.GridExp.filename = filenamegrid;
    handles.GridExp.zone = zone;
end
handles = PlotMultiSessionIn1(handles);
if handles.VortLate.V ==1
    SideVortS(handles)
end
display('All done!')


% --- Executes on selection change in VelPopUp.
function VelPopUp_Callback(hObject, eventdata, handles)
% hObject    handle to VelPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = get(hObject,'Value');
GraficosPPlotar = handles.GraficosPPlotar;
handles.VelPopUp.V = GraficosPPlotar{contents};
guidata(hObject,handles)
% Hints: contents = cellstr(get(hObject,'String')) returns VelPopUp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from VelPopUp


% --- Executes during object creation, after setting all properties.
function VelPopUp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VelPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Chen = cheneq(varargin)
source = varargin{1};
Chen = get(source,'Label');
parent = get(source,'Parent');
data = get(parent,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
switch Chen
    case 'Z_0'
        input = inputdlg('Value for Chens equation','Input',1,datastr(2));
        if sum(size(input)==0)~=0
            return
        end
        
        data(2) = str2num(input{1});
        set(parent,'UserData',data);
    case 'U*'
        input = inputdlg('Value for Chens equation','Input',1,datastr(1));
        if sum(size(input)==0)~=0
            return
        end
        data(1) = str2num(input{1});
        set(parent,'UserData',data);   
    case 'n'
        input = inputdlg('Value for Chens equation','Input',1,datastr(3));
        if sum(size(input)==0)~=0
            return
        end
        data(3) = str2num(input{1});
        set(parent,'UserData',data);    
end

function numbat = Inputnumbat(varargin)
data = get(varargin{1},'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
numbat = inputdlg('Number of points/scale','Input',1,datastr(1));
if sum(size(numbat)==0)~=0
    return
end
numbat = str2num(numbat{1});
set(varargin{1},'UserData',numbat);

function shp2d = InputExpShp2d(varargin)
source = varargin{1};
shp2d = get(source,'Label');
parent = get(source,'Parent');
switch shp2d
    case 'Geopoint'
        set(parent,'UserData','Geopoint');
    case 'Geovector'
        set(parent,'UserData','Geovector');        
end

function shpfileexp = shpfileexpfun(varargin)
source = varargin{1};
shpfileexp = get(source,'Label');
parent = get(source,'Parent');
switch shpfileexp
    case 'utm'
        set(parent,'UserData','utm');
    case 'Lat & Long'
        set(parent,'UserData','Lat & Long');        
end

function plotpos = PlotPos(varargin)
data = get(varargin{1},'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
plp = inputdlg({'Distance from left edge of monitor','Distance from bottom edge of monitor',...
    'Length','Height'},'Input',1,{datastr{1:4}});
if sum(size(plp)==0)~=0
    return
end
plp1 = str2num(plp{1});
plp2 = str2num(plp{2});
plp3 = str2num(plp{3});
plp4 = str2num(plp{4});
plotpos = [plp1 plp2 plp3 plp4];
set(varargin{1},'UserData',plotpos);

function sideVErtout = SideVert(varargin)
source = varargin{1};
sideVErtout = get(source,'Label');
parent = get(source,'Parent');
data = get(parent,'UserData');
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
switch sideVErtout
    case 'Number of Arrows'
        input = inputdlg('Number of Arrows:','Input',1,datastr(1));
        if sum(size(input)==0)~=0
            return
        end
        input = str2num(input{1});
        set(parent,'UserData',[input data(2)]);
    case 'Percentage of the length'
        input = inputdlg('Percentage of the Length:','Input',1,datastr(2));
        if sum(size(input)==0)~=0
            return
        end
        input = str2num(input{1});
        if or(input<0,input>1)
            error('Percentage must be bigger than 0 and smaller than 1')
        end
        set(parent,'UserData',[data(1) input]);        
end
