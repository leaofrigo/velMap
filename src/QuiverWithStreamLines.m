function varargout = QuiverWithStreamLines(varargin)
% QUIVERWITHSTREAMLINES MATLAB code for QuiverWithStreamLines.fig
%      QUIVERWITHSTREAMLINES, by itself, creates a new QUIVERWITHSTREAMLINES or raises the existing
%      singleton*.
%
%      H = QUIVERWITHSTREAMLINES returns the handle to a new QUIVERWITHSTREAMLINES or the handle to
%      the existing singleton*.
%
%      QUIVERWITHSTREAMLINES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in QUIVERWITHSTREAMLINES.M with the given input arguments.
%
%      QUIVERWITHSTREAMLINES('Property','Value',...) creates a new QUIVERWITHSTREAMLINES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before QuiverWithStreamLines_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to QuiverWithStreamLines_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help QuiverWithStreamLines

% Last Modified by GUIDE v2.5 02-Jun-2015 16:26:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @QuiverWithStreamLines_OpeningFcn, ...
                   'gui_OutputFcn',  @QuiverWithStreamLines_OutputFcn, ...
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


% --- Executes just before QuiverWithStreamLines is made visible.
function QuiverWithStreamLines_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to QuiverWithStreamLines (see VARARGIN)
handles.properties.Arrowsize = 1;
handles.properties.colour = 'b';
handles.properties.colourstream = 'b';
handles.properties.line = .5;
n = handles.properties.Arrowsize;
colour = handles.properties.colour;
GridVar = varargin{1};
numofvalues = 100;
x = linspace(min(GridVar.Longitudeidx),max(GridVar.Longitudeidx),numofvalues);
y = linspace(min(GridVar.Latitudeidx),max(GridVar.Latitudeidx),numofvalues);
[X,Y] = meshgrid(x,y);
Ve = griddata(GridVar.Longitudeidx,GridVar.Latitudeidx,GridVar.VelocityEidx,X,Y);
Vn = griddata(GridVar.Longitudeidx,GridVar.Latitudeidx,GridVar.VelocityNidx,X,Y);
% quiver(X,Y,Ve,Vn)
QuivHandle = quiver(GridVar.Longitudeidx,GridVar.Latitudeidx,GridVar.VelocityEidx,GridVar.VelocityNidx,...
    n,'Color',colour);
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
xlim([min(X(:)) max(X(:))]);
ylim([min(Y(:)) max(Y(:))]);
handles.GridVar=GridVar;
handles.x=x;
handles.y=y;
handles.X=X;
handles.Y=Y;
handles.Ve=Ve;
handles.Vn=Vn;
data = [];
handles.QuivHandle = QuivHandle;
handles.data= data;
set(handles.uitable1,'data',data,'ColumnName',{'X','Y','S/N'},'ColumnFormat',{'numeric','numeric','logical'},...
    'ColumnEditable', [false false true],'RowName',[],'ColumnWidth',{80, 80,80});
% Choose default command line output for QuiverWithStreamLines
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes QuiverWithStreamLines wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = QuiverWithStreamLines_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = handles.data;
[posx,posy] = ginput;
datatemp=[];
for i=1:length(posx)
    datatemp = [datatemp;{posx(i) posy(i) true}];
end
data = [data;datatemp];
handles.data = data;
set(handles.uitable1,'data',data)
guidata(hObject, handles);


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.data = [];
data = handles.data;
set(handles.uitable1,'data',data);
guidata(hObject,handles)


% --- Executes when selected cell(s) is changed in uitable1.
function uitable1_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
% editData = eventdata.EditData;
% indData = eventdata.Indices;
% data = handles.data;
% tabledata{indData(1),2} = editData;
% handles.data = data;
% set(handles.uitable1,'data',handles.data)
% guidata(hObject,handles)


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
data= handles.data;
% if editData == 0
%     editData = 'false';
% elseif editData==1
%     editData = 'true';
% end
data{indData(1),indData(2)} = editData;
handles.data = data;
set(handles.uitable1,'data',data);
guidata(hObject,handles)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
X = handles.X;
Y = handles.Y;
Ve= handles.Ve;
Vn= handles.Vn;
data = handles.data;
startx = [data{:,1}];
starty = [data{:,2}];
tf = [data{:,3}];
startx = startx(tf); starty = starty(tf);
h = streamline(X,Y,Ve,Vn,startx,starty);
% save('C:\Users\Roberta\Desktop\Ricardo\Variables Save\handles.mat','handles')
for i =1:length(h)
    set(h,'Color',handles.properties.colourstream)
    set(h,'LineWidth',handles.properties.line)
end
handles.handleStream = h;
guidata(hObject,handles)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton4.
function pushbutton4_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in clearplot.
function clearplot_Callback(hObject, eventdata, handles)
% hObject    handle to clearplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
X = handles.X;
Y = handles.Y;
Ve = handles.Ve;
Vn = handles.Vn;
cla(handles.axes1)
GridVar = handles.GridVar;
n = handles.properties.Arrowsize;
colour = handles.properties.colour;
% save('C:\Users\Roberta\Desktop\Ricardo\Variables Save\handles.mat','handles')

handles = rmfield(handles,'QuivHandle');
A=any(strcmp('handleStream',fieldnames(handles)));
if A~=0
    handles = rmfield(handles,'handleStream');
end
QuivHandle = quiver(GridVar.Longitudeidx,GridVar.Latitudeidx,GridVar.VelocityEidx,GridVar.VelocityNidx,...
    n,'Color',colour);
handles.QuivHandle = QuivHandle;
xlim([min(X(:)) max(X(:))]);
ylim([min(Y(:)) max(Y(:))]);
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

guidata(hObject,handles)


% --- Executes on button press in Graphproperties.
function Graphproperties_Callback(hObject, eventdata, handles)
% hObject    handle to Graphproperties (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% properties = PlotPropertiesStreamLines;
% handles.properties = properties;
% clearplot_Callback(hObject, eventdata, handles)
% guidata(hObject,handles)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over Graphproperties.
function Graphproperties_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Graphproperties (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over arrowsize.
function arrowsize_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to arrowsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = handles.properties.Arrowsize;
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
prompt = {'Input Arrow Size Multiplier (default = 1)'};
Arrowsize = inputdlg(prompt,'Input',1,datastr);
if sum(size(Arrowsize)==0)~=0
    return
end
Arrowsize = str2num(Arrowsize{1});
handles.properties.Arrowsize = Arrowsize;
h = handles.QuivHandle;
for i =1:length(h)
    set(h,'AutoScaleFactor',handles.properties.Arrowsize)
end
guidata(hObject,handles)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over quivcolor.
function quivcolor_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to quivcolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prompt = {'Choose a Color (default = blue):'};
ColorsAvailable = {'yellow'; 'magenta'; 'cyan';'green';'blue';'white';'black'};
[colour, ok]= listdlg('PromptString',prompt,'SelectionMode','single','ListString',ColorsAvailable,...
    'ListSize',[160 105]);
if ok ==0
    return
end
if colour == 1
    colour = 'y';
elseif colour == 2
    colour = 'm';
elseif colour == 3
    colour = 'c';
elseif colour == 4
    colour = 'g';
elseif colour == 5
    colour = 'b';
elseif colour ==6
    colour = 'w';
elseif colour ==7
    colour = 'k';
end
handles.properties.colour = colour;
h = handles.QuivHandle;
for i =1:length(h)
    set(h,'Color',handles.properties.colour)
end
guidata(hObject,handles)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over ColorStream.
function ColorStream_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to ColorStream (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prompt = {'Choose a Color (default = blue):'};
ColorsAvailable = {'yellow'; 'magenta'; 'cyan';'green';'blue';'white';'black'};
colour = listdlg('PromptString',prompt,'SelectionMode','single','ListString',ColorsAvailable,...
    'ListSize',[160 105]);
% colour = ColorsAvailable(v);
if colour == 1
    colour = 'y';
elseif colour == 2
    colour = 'm';
elseif colour == 3
    colour = 'c';
elseif colour == 4
    colour = 'g';
elseif colour == 5
    colour = 'b';
elseif colour ==6
    colour = 'w';
elseif colour ==7
    colour = 'k';
end
handles.properties.colourstream = colour;
A=any(strcmp('handleStream',fieldnames(handles)));
if A == 0
    guidata(hObject,handles)
    return
end
h=handles.handleStream;
% handles.properties.colourstream = colour;
% A = exist('handles.handleStream','var');
% if A == 0
%     guidata(hObject,handles)
%     return
% end
for i =1:length(h)
    set(h,'Color',handles.properties.colourstream)
    set(h,'LineWidth',handles.properties.line)
end
guidata(hObject,handles)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over LineWidthStream.
function LineWidthStream_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to LineWidthStream (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = handles.properties.line;
datastr = cellfun(@num2str, num2cell(data), 'UniformOutput', false);
prompt = {'Input StreamLine width (default = 0.5)'};
line = inputdlg(prompt,'Input',1,datastr(1));
if sum(size(line)==0)~=0
    return
end
line = line{1};
handles.properties.line = str2num(line);
A=any(strcmp('handleStream',fieldnames(handles)));
if A == 0
    guidata(hObject,handles)
    return
end
h = handles.handleStream;
for i =1:length(h)
    set(h,'LineWidth',handles.properties.line)
end
guidata(hObject,handles)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over reset.
function reset_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.properties.Arrowsize = 1;
handles.properties.colour = 'b';
handles.properties.colourstream = 'b';
handles.properties.line = .5;
A = exist('handles.handleStream','var');
if A == 0
    h1 = handles.QuivHandle;
    for i =1:length(h1)
        set(h1,'Color',handles.properties.colour)
        set(h1,'AutoScaleFactor',handles.properties.Arrowsize)
    end    
    guidata(hObject,handles)
    return
end
h = handles.handleStream;
for i =1:length(h)
    set(h,'LineWidth',handles.properties.line)
    set(h,'Color',handles.properties.colourstream)
end
h1 = handles.QuivHandle;
for i =1:length(h1)
    set(h1,'Color',handles.properties.colour)
    set(h1,'AutoScaleFactor',handles.properties.Arrowsize)
end
guidata(hObject,handles)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over copyfigure.
function copyfigure_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to copyfigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h1 = figure('Visible','on');
X = handles.X;
Y = handles.Y;
Ve = handles.Ve;
Vn = handles.Vn;
GridVar = handles.GridVar;
n = handles.properties.Arrowsize;
colour = handles.properties.colour;
quiver(GridVar.Longitudeidx,GridVar.Latitudeidx,GridVar.VelocityEidx,GridVar.VelocityNidx,...
    n,'Color',colour);
data = handles.data;
startx = [data{:,1}];
starty = [data{:,2}];
tf = [data{:,3}];
startx = startx(tf); starty = starty(tf);
h = streamline(X,Y,Ve,Vn,startx,starty);
for i =1:length(h)
    set(h,'LineWidth',handles.properties.line)
    set(h,'Color',handles.properties.colourstream)
end
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
xlim([min(X(:)) max(X(:))]);
ylim([min(Y(:)) max(Y(:))]);
% Pos = get(handles.axes1,'OuterPosition');
% get(handles.axes1)
% set(h1,'Position',Pos)

%plot_google_map
print(h1,'-dmeta')
savefig(h1,[pwd '\StreamLineQuiver.fig'])
%close(h1)


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
