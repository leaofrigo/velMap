function varargout = quiver2DGUI(varargin)
% QUIVER2DGUI MATLAB code for quiver2DGUI.fig
%      QUIVER2DGUI, by itself, creates a new QUIVER2DGUI or raises the existing
%      singleton*.
%
%      H = QUIVER2DGUI returns the handle to a new QUIVER2DGUI or the handle to
%      the existing singleton*.
%
%      QUIVER2DGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in QUIVER2DGUI.M with the given input arguments.
%
%      QUIVER2DGUI('Property','Value',...) creates a new QUIVER2DGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before quiver2DGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to quiver2DGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help quiver2DGUI

% Last Modified by GUIDE v2.5 18-Jun-2015 16:09:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @quiver2DGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @quiver2DGUI_OutputFcn, ...
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


% --- Executes just before quiver2DGUI is made visible.
function quiver2DGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to quiver2DGUI (see VARARGIN)
Velocidades = [{'North'} {'East'} {'Upstream'} {'Magnitude'} {'Error'} {'Transverse'} {'Longitudinal'}];
set(handles.Velocidade1,'String',Velocidades)
set(handles.Velocidade2,'String',Velocidades)
handles.Vel1V = 'North';
handles.Vel2V = 'North';
% Choose default command line output for quiver2DGUI
% handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes quiver2DGUI wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = quiver2DGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Vel1 = handles.Vel1V;
Vel2 = handles.Vel2V;
varargout = [{Vel1} {Vel2}];
handles.output = varargout;
delete(handles.figure1)
% Get default command line output from handles structure
% varargout{1} = handles.output;


% --- Executes on selection change in Velocidade1.
function Velocidade1_Callback(hObject, eventdata, handles)
% hObject    handle to Velocidade1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = get(hObject,'Value');
Vels =[{'North'} {'East'} {'Upstream'} {'Magnitude'} {'Error'} {'Transverse'} {'Longitudinal'}];
handles.Vel1V = Vels{contents};
guidata(hObject,handles)
% Hints: contents = cellstr(get(hObject,'String')) returns Velocidade1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Velocidade1


% --- Executes during object creation, after setting all properties.
function Velocidade1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Velocidade1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Velocidade2.
function Velocidade2_Callback(hObject, eventdata, handles)
% hObject    handle to Velocidade2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = get(hObject,'Value');
Vels =[{'North'} {'East'} {'Upstream'} {'Magnitude'} {'Error'} {'Transverse'} {'Longitudinal'}];
handles.Vel2V = Vels{contents};
guidata(hObject,handles)
% Hints: contents = cellstr(get(hObject,'String')) returns Velocidade2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Velocidade2


% --- Executes during object creation, after setting all properties.
function Velocidade2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Velocidade2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Pronto.
function Pronto_Callback(hObject, eventdata, handles)
% hObject    handle to Pronto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1)
% quiver2DGUI_OutputFcn(hObject, eventdata, handles)
