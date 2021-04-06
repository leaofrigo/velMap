function varargout = SimuVelocidades(varargin)
% SIMUVELOCIDADES MATLAB code for SimuVelocidades.fig
%      SIMUVELOCIDADES, by itself, creates a new SIMUVELOCIDADES or raises the existing
%      singleton*.
%
%      H = SIMUVELOCIDADES returns the handle to a new SIMUVELOCIDADES or the handle to
%      the existing singleton*.
%
%      SIMUVELOCIDADES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SIMUVELOCIDADES.M with the given input arguments.
%
%      SIMUVELOCIDADES('Property','Value',...) creates a new SIMUVELOCIDADES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SimuVelocidades_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SimuVelocidades_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SimuVelocidades

% Last Modified by GUIDE v2.5 19-May-2015 15:51:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SimuVelocidades_OpeningFcn, ...
                   'gui_OutputFcn',  @SimuVelocidades_OutputFcn, ...
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


% --- Executes just before SimuVelocidades is made visible.
function SimuVelocidades_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SimuVelocidades (see VARARGIN)
Velocidades = [{'North'} {'East'} {'Upstream'} {'Magnitude'}];
set(handles.Velocidades,'String',Velocidades)
handles.Vel1V = 'North';
% Choose default command line output for SimuVelocidades
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SimuVelocidades wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SimuVelocidades_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
Vel1 = handles.Vel1V;
varargout = [{Vel1}];
handles.output = varargout;
delete(handles.figure1)


% --- Executes on selection change in Velocidades.
function Velocidades_Callback(hObject, eventdata, handles)
% hObject    handle to Velocidades (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = get(hObject,'Value');
Vels =[{'North'} {'East'} {'Upstream'} {'Magnitude'}];
handles.Vel1V = Vels{contents};
guidata(hObject,handles)
% Hints: contents = cellstr(get(hObject,'String')) returns Velocidades contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Velocidades


% --- Executes during object creation, after setting all properties.
function Velocidades_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Velocidades (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1)
