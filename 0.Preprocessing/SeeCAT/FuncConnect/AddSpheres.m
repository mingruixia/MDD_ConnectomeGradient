function varargout = AddSpheres(varargin)
% ADDSPHERES MATLAB code for AddSpheres.fig
%      ADDSPHERES, by itself, creates a new ADDSPHERES or raises the existing
%      singleton*.
%
%      H = ADDSPHERES returns the handle to a new ADDSPHERES or the handle to
%      the existing singleton*.
%
%      ADDSPHERES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ADDSPHERES.M with the given input arguments.
%
%      ADDSPHERES('Property','Value',...) creates a new ADDSPHERES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AddSpheres_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AddSpheres_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AddSpheres

% Last Modified by GUIDE v2.5 02-Jul-2016 14:57:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AddSpheres_OpeningFcn, ...
                   'gui_OutputFcn',  @AddSpheres_OutputFcn, ...
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


% --- Executes just before AddSpheres is made visible.
function AddSpheres_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AddSpheres (see VARARGIN)

% Choose default command line output for AddSpheres
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

movegui(hObject,'center');

% UIWAIT makes AddSpheres wait for user response (see UIRESUME)
uiwait(handles.AddSpheres_figure);



% --- Outputs from this function are returned to the command line.
function varargout = AddSpheres_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if isempty(handles)
    varargout{1}=[];
else
    varargout{1} = handles.CoordMat;
    delete(handles.AddSpheres_figure);
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in T2M_checkbox.
function T2M_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to T2M_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of T2M_checkbox


% --- Executes on button press in OK_pushbutton.
function OK_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to OK_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CoordMat = str2num(get(handles.Seed_edit, 'String'));
if get(handles.T2M_checkbox, 'Value') == 1
    CoordMat = [seecat_tal2icbm_spm(CoordMat(:,1:3)),CoordMat(:,4)];
end
handles.CoordMat = CoordMat;
guidata(hObject,handles);
uiresume(handles.AddSpheres_figure);

function outpoints = seecat_tal2icbm_spm(inpoints)
%
% This function converts coordinates from Talairach space to MNI
% space (normalized using the SPM software package) using the 
% tal2icbm transform developed and validated by Jack Lancaster 
% at the Research Imaging Center in San Antonio, Texas.
%
% http://www3.interscience.wiley.com/cgi-bin/abstract/114104479/ABSTRACT
% 
% FORMAT outpoints = icbm_spm2tal(inpoints)
% Where inpoints is N by 3 or 3 by N matrix of coordinates
% (N being the number of points)
%
% ric.uthscsa.edu 3/14/07

% find which dimensions are of size 3
dimdim = find(size(inpoints) == 3);
if isempty(dimdim)
  error('input must be a N by 3 or 3 by N matrix')
end

% 3x3 matrices are ambiguous
% default to coordinates within a row
if dimdim == [1 2]
  disp('input is an ambiguous 3 by 3 matrix')
  disp('assuming coordinates are row vectors')
  dimdim = 2;
end

% transpose if necessary
if dimdim == 2
  inpoints = inpoints';
end

% Transformation matrices, different for each software package
icbm_spm = [0.9254 0.0024 -0.0118 -1.0207
	   	   -0.0048 0.9316 -0.0871 -1.7667
            0.0152 0.0883  0.8924  4.0926
            0.0000 0.0000  0.0000  1.0000];

% invert the transformation matrix
icbm_spm = inv(icbm_spm);

% apply the transformation matrix
inpoints = [inpoints; ones(1, size(inpoints, 2))];
inpoints = icbm_spm * inpoints;

% format the outpoints, transpose if necessary
outpoints = inpoints(1:3, :);
if dimdim == 2
  outpoints = outpoints';
end




% --- Executes on button press in Cancel_pushbutton.
function Cancel_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.AddSpheres_figure);



function Seed_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Seed_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Seed_edit as text
%        str2double(get(hObject,'String')) returns contents of Seed_edit as a double


% --- Executes during object creation, after setting all properties.
function Seed_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Seed_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
