function varargout = GUI_Trace(varargin)
% GUI_TRACE MATLAB code for GUI_Trace.fig
%      GUI_TRACE, by itself, creates a new GUI_TRACE or raises the existing
%      singleton*.
%
%      H = GUI_TRACE returns the handle to a new GUI_TRACE or the handle to
%      the existing singleton*.
%
%      GUI_TRACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_TRACE.M with the given input arguments.
%
%      GUI_TRACE('Property','Value',...) creates a new GUI_TRACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_Trace_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_Trace_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_Trace

% Last Modified by GUIDE v2.5 15-Aug-2018 14:22:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;

gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Trace_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Trace_OutputFcn, ...
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





% --- Executes just before GUI_Trace is made visible.
function GUI_Trace_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_Trace (see VARARGIN)

% Choose default command line output for GUI_Trace
handles.output = hObject;


% Update handles structure
guidata(hObject, handles);
global settings_plt;
settings_plt{1,1} = hObject;

% UIWAIT makes GUI_Trace wait for user response (see UIRESUME)
% uiwait(handles.Settings);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_Trace_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
waitfor(hObject,'Tag');
varargout{1} = str2double(get(handles.brush,'String'));
varargout{2} = str2double(get(handles.width,'String'));
varargout{3} = str2double(get(handles.height,'String'));
varargout{4} = get(handles.radiobutton_pp,'Value');
varargout{5} = str2double(get(handles.ncolors,'String'));
close(hObject);



% --- Executes on button press in ok.
function ok_Callback(hObject, eventdata, handles)
% hObject    handle to ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[h,figure] = gcbo;
figure.Tag = "Saving...";



% --- Executes on button press in radiobutton_pp.
function radiobutton_pp_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_pp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_pp
global settings_plt;
settings_plt{4,1} = 1;


% --- Executes on button press in radiobutton_pc.
function radiobutton_pc_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_pc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_pc
global settings_plt;
settings_plt{4,1} = 0;



function width_Callback(hObject, eventdata, handles)
% hObject    handle to width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of width as text
%        str2double(get(hObject,'String')) returns contents of width as a double
global settings_plt;
settings_plt{2,1} = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function height_Callback(hObject, eventdata, handles)
% hObject    handle to height (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of height as text
%        str2double(get(hObject,'String')) returns contents of height as a double
global settings_plt;
settings_plt{3,1} = str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function height_CreateFcn(hObject, eventdata, handles)
% hObject    handle to height (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close Settings.
function Settings_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to Settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);



function brush_Callback(hObject, eventdata, handles)
% hObject    handle to brush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of brush as text
%        str2double(get(hObject,'String')) returns contents of brush as a double


% --- Executes during object creation, after setting all properties.
function brush_CreateFcn(hObject, eventdata, handles)
% hObject    handle to brush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ncolors_Callback(hObject, eventdata, handles)
% hObject    handle to ncolors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ncolors as text
%        str2double(get(hObject,'String')) returns contents of ncolors as a double


% --- Executes during object creation, after setting all properties.
function ncolors_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ncolors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
