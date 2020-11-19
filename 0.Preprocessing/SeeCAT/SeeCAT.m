function varargout = SeeCAT(varargin)
% SeeCAT (Seed-based Connectivity Analysis Toolbox) is designed for
% resting-state fMRI data preprocessing (with SPM), seed based functional
% and structural connectivity computation, voxel-based functional
% connectivity strength calculation, and statistic and visulization for
% volume based metrics.
%-----------------------------------------------------------
%	Copyright(c) 2017
%	Beijing Normal University
%	Written by Mingrui Xia (Preprocessing for fMRI data, Seed-based 
%   connectivity, and voxel-based connectivity strength) and Xindi Wang 
%   (Statistical analysis and viewer)
%	Mail to Author: <a href="mingruixia@gmail.com">Mingrui Xia</a> or <a href="sandywang.rest@live.com">Xindi Wang</a>
%   Version 1.0 alpha;
%   First generated 201604012;
%   Last edited 20170705
%-----------------------------------------------------------
%
% SEECAT MATLAB code for SeeCAT.fig
%      SEECAT, by itself, creates a new SEECAT or raises the existing
%      singleton*.
%
%      H = SEECAT returns the handle to a new SEECAT or the handle to
%      the existing singleton*.
%
%      SEECAT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEECAT.M with the given input arguments.
%
%      SEECAT('Property','Value',...) creates a new SEECAT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SeeCAT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SeeCAT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SeeCAT

% Last Modified by GUIDE v2.5 01-Jul-2016 22:54:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SeeCAT_OpeningFcn, ...
                   'gui_OutputFcn',  @SeeCAT_OutputFcn, ...
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


% --- Executes just before SeeCAT is made visible.
function SeeCAT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SeeCAT (see VARARGIN)

% Choose default command line output for SeeCAT
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SeeCAT wait for user response (see UIRESUME)
% uiwait(handles.SeeCAT_figure);
movegui(hObject,'center');

% --- Outputs from this function are returned to the command line.
function varargout = SeeCAT_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in PrefMRI_pushbutton.
function PrefMRI_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to PrefMRI_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SeeCAT_PrepfMRI;


% --- Executes on button press in FuncConnect_pushbutton.
function FuncConnect_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to FuncConnect_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SeeCAT_FuncConnect;


% --- Executes on button press in VoxDeg_pushbutton.
function VoxDeg_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to VoxDeg_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SeeCAT_VoxDeg;


% --- Executes on button press in Stat_pushbutton.
function Stat_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Stat_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SeeCAT_Stat;


% --- Executes on button press in Visual_pushbutton.
function Visual_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Visual_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SeeCAT_Viewer;


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ImgCal_pushbutton.
function ImgCal_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to ImgCal_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SeeCAT_Calculator;
