function varargout = PlotPropertiesStreamLines(varargin)
% PLOTPROPERTIESSTREAMLINES MATLAB code for PlotPropertiesStreamLines.fig
%      PLOTPROPERTIESSTREAMLINES, by itself, creates a new PLOTPROPERTIESSTREAMLINES or raises the existing
%      singleton*.
%
%      H = PLOTPROPERTIESSTREAMLINES returns the handle to a new PLOTPROPERTIESSTREAMLINES or the handle to
%      the existing singleton*.
%
%      PLOTPROPERTIESSTREAMLINES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOTPROPERTIESSTREAMLINES.M with the given input arguments.
%
%      PLOTPROPERTIESSTREAMLINES('Property','Value',...) creates a new PLOTPROPERTIESSTREAMLINES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PlotPropertiesStreamLines_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PlotPropertiesStreamLines_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PlotPropertiesStreamLines

% Last Modified by GUIDE v2.5 22-May-2015 13:33:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PlotPropertiesStreamLines_OpeningFcn, ...
                   'gui_OutputFcn',  @PlotPropertiesStreamLines_OutputFcn, ...
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


% --- Executes just before PlotPropertiesStreamLines is made visible.
function PlotPropertiesStreamLines_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PlotPropertiesStreamLines (see VARARGIN)

% Choose default command line output for PlotPropertiesStreamLines
handles.output = hObject;
handles.colourStream = 'b';
handles.colour = 'b';
handles.Arrowsize = '1';
handles.line = '0.5';
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PlotPropertiesStreamLines wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PlotPropertiesStreamLines_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Prop.ArrowSize = handles.Arrowsize;
Prop.colour = handles.colour;
Prop.colourstream = handles.colourStream;
Prop.line = handles.line;
varargout = [{Prop}];
handles.output = varargout;
delete(handles.figure1)
% Get default command line output from handles structure


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text1.
function text1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to text1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prompt = {'Input Arrow Size Multiplier (default = 1)'};
Arrowsize = inputdlg(prompt);
Arrowsize = str2mat(Arrowsize);
handles.Arrowsize = Arrowsize;
guidata(hObject,handles)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over Colour.
function Colour_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Colour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prompt = {'Choose a Color (default = blue):'};
ColorsAvailable = {'yellow'; 'magenta'; 'cyan';'green';'blue';'white';'black'};
colour = listdlg('PromptString',prompt,'SelectionMode','single','ListString',ColorsAvailable,...
    'ListSize',[160 105]);
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
handles.colour = colour;
guidata(hObject,handles)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over ColourStream.
function ColourStream_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to ColourStream (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prompt = {'Choose a Color (default = blue):'};
ColorsAvailable = {'yellow'; 'magenta'; 'cyan';'green';'blue';'white';'black'};
colour = listdlg('PromptString',prompt,'SelectionMode','single','ListString',ColorsAvailable,...
    'ListSize',[160 105]);
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
handles.colourStream = colour;
guidata(hObject,handles)

% --- Executes on button press in Close.
function Close_Callback(hObject, eventdata, handles)
% hObject    handle to Close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text4.
function text4_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to text4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prompt = {'Input StreamLine width (default = 0.5)'};
line = inputdlg(prompt);
line = line{1};
handles.line = line;
guidata(hObject,handles)
