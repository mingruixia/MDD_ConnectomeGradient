function varargout = SeeCAT_Viewer(varargin)
% SeeCAT_Viewer M-file for SeeCAT_Viewer.fig
%      SeeCAT_Viewer, by itself, creates a new SeeCAT_Viewer or raises the existing
%      singleton*.
%
%      H = SeeCAT_Viewer returns the handle to a new SeeCAT_Viewer or the handle to
%      the existing singleton*.
%
%      SeeCAT_Viewer('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SeeCAT_Viewer.M with the given input arguments.
%
%      SeeCAT_Viewer('Property','Value',...) creates a new SeeCAT_Viewer or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SeeCAT_Viewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SeeCAT_Viewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Written by Wang Xin-di, 20130809
% State Key Laboratory of Cognitive Neuroscience and Learning, Beijing
% Normal University, Beijing, PR China
% sandywang.rest@gmail.com

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SeeCAT_Viewer_OpeningFcn, ...
                   'gui_OutputFcn',  @SeeCAT_Viewer_OutputFcn, ...
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


% --- Executes just before SeeCAT_Viewer is made visible.
function SeeCAT_Viewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SeeCAT_Viewer (see VARARGIN)

% Choose default command line output for SeeCAT_Viewer
movegui(hObject, 'center');

addlistener(handles.ThrdSlider, 'Value',...
    'PostSet', ...
    @(objH, eventData) ThrdListener_Callback(objH, eventData, hObject));
addlistener(handles.AlphaSlider, 'Value',...
    'PostSet', ...
    @(objH, eventData) AlphaListener_Callback(objH, eventData, hObject));


[DPABIPath, fileN, extn] = fileparts(which('SeeCAT.m'));
TemplatePath=fullfile(DPABIPath, 'templates');

NiiFile=fullfile(TemplatePath, 'aal.nii');
MatFile=fullfile(TemplatePath, 'aal_Labels.mat');
AALInfo=w_GetAtlasInfo(MatFile, NiiFile, 'AAL');
NiiFile=fullfile(TemplatePath, 'Brodmann_YCG.nii');
MatFile=fullfile(TemplatePath, 'Brodmann_YCG_Labels.mat');
BrodInfo=w_GetAtlasInfo(MatFile, NiiFile, 'Brodmann');
AtlasInfo=cell(2, 1);
AtlasInfo{1, 1}=AALInfo;
AtlasInfo{2, 1}=BrodInfo;

global st

if ~iscell(st)
    st=cell(24);
end
curfig=w_Compatible2014bFig(hObject);

st{curfig}.fig=hObject;

st{curfig}.curblob=0;
st{curfig}.AtlasInfo=AtlasInfo;
st{curfig}.xhairs=1;st{curfig}.hld=1;st{curfig}.yoke=0;
st{curfig}.n=0;st{curfig}.vols=cell(2);st{curfig}.bb=[];
st{curfig}.Space=eye(4,4);st{curfig}.centre=[0 0 0];
st{curfig}.callback=';';st{curfig}.mode=1;st{curfig}.snap=[];
st{curfig}.plugins={};%'reorient' 'roi'   'rgb'  
st{curfig}.TCFlag=[];st{curfig}.SSFlag=[];
st{curfig}.MPFlag=[];

colormap(gray(64));

UnderlayFileName=fullfile(TemplatePath,'ch2.nii');

[UnderlayVolume UnderlayVox UnderlayHeader] = y_ReadRPI(UnderlayFileName);
handles.UnderlayFileName='';
handles.UserDefinedFileName=UnderlayFileName;
set(handles.UnderlayEntry, 'String', 'ch2.nii');

UnderlayHeader.Data = UnderlayVolume;
UnderlayHeader.Vox  = UnderlayVox;
%UnderlayHeader.rmat = UnderlayHeader.mat;

handles.output = hObject;
handles.OverlayHeaders=cell(5, 1);

% Update handles structure
guidata(hObject, handles);

w_spm_orthviews('Image',UnderlayHeader);
%w_spm_orthviews('AddContext',1);

% UIWAIT makes SeeCAT_Viewer wait for user response (see UIRESUME)
% uiwait(handles.DPABI_fig);


% --- Outputs from this function are returned to the command line.
function varargout = SeeCAT_Viewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function XEntry_Callback(hObject, eventdata, handles)
% hObject    handle to XEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
centre=Pos(handles, 'XYZ');
w_spm_orthviews('Reposition', centre);
% Hints: get(hObject,'String') returns contents of XEntry as text
%        str2double(get(hObject,'String')) returns contents of XEntry as a double


% --- Executes during object creation, after setting all properties.
function XEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function YEntry_Callback(hObject, eventdata, handles)
% hObject    handle to YEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
centre=Pos(handles, 'XYZ');
w_spm_orthviews('Reposition', centre);
% Hints: get(hObject,'String') returns contents of YEntry as text
%        str2double(get(hObject,'String')) returns contents of YEntry as a double


% --- Executes during object creation, after setting all properties.
function YEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to YEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ZEntry_Callback(hObject, eventdata, handles)
% hObject    handle to ZEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
centre=Pos(handles, 'XYZ');
w_spm_orthviews('Reposition', centre);
% Hints: get(hObject,'String') returns contents of ZEntry as text
%        str2double(get(hObject,'String')) returns contents of ZEntry as a double


% --- Executes during object creation, after setting all properties.
function ZEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ZEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function KEntry_Callback(hObject, eventdata, handles)
% hObject    handle to KEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
centre=Pos(handles, 'IJK');
w_spm_orthviews('Reposition', centre);
% Hints: get(hObject,'String') returns contents of KEntry as text
%        str2double(get(hObject,'String')) returns contents of KEntry as a double


% --- Executes during object creation, after setting all properties.
function KEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to KEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function JEntry_Callback(hObject, eventdata, handles)
% hObject    handle to JEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
centre=Pos(handles, 'IJK');
w_spm_orthviews('Reposition', centre);
% Hints: get(hObject,'String') returns contents of JEntry as text
%        str2double(get(hObject,'String')) returns contents of JEntry as a double


% --- Executes during object creation, after setting all properties.
function JEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to JEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function IEntry_Callback(hObject, eventdata, handles)
% hObject    handle to IEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
centre=Pos(handles, 'IJK');
w_spm_orthviews('Reposition', centre);
% Hints: get(hObject,'String') returns contents of IEntry as text
%        str2double(get(hObject,'String')) returns contents of IEntry as a
%        double

% --- Executes during object creation, after setting all properties.
function IEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in AButton.
function AButton_Callback(hObject, eventdata, handles)
% hObject    handle to AButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global st
curfig=handles.DPABI_fig;
curfig=w_Compatible2014bFig(curfig);

TMFlag=w_Montage(handles.DPABI_fig, 'T');
%TMFlag=w_Compatible2014bFig(TMFlag);
st{curfig}.MPFlag=[st{curfig}.MPFlag;{TMFlag}];

% --- Executes on button press in SButton.
function SButton_Callback(hObject, eventdata, handles)
% hObject    handle to SButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global st
curfig=handles.DPABI_fig;
curfig=w_Compatible2014bFig(curfig);

SMFlag=w_Montage(handles.DPABI_fig, 'S');
%SMFlag=w_Compatible2014bFig(SMFlag);
st{curfig}.MPFlag=[st{curfig}.MPFlag;{SMFlag}];

% --- Executes on button press in CButton.
function CButton_Callback(hObject, eventdata, handles)
% hObject    handle to CButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global st
curfig=handles.DPABI_fig;
curfig=w_Compatible2014bFig(curfig);

CMFlag=w_Montage(handles.DPABI_fig, 'C');
%CMFlag=w_Compatible2014bFig(CMFlag);
st{curfig}.MPFlag=[st{curfig}.MPFlag;{CMFlag}];

function ScaleEntry_Callback(hObject, eventdata, handles)
% hObject    handle to ScaleEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global st
curfig=handles.DPABI_fig;
curfig=w_Compatible2014bFig(curfig);

scale=get(handles.ScaleEntry, 'String');
scale=str2double(scale);

if isempty(scale)
    scale=100;
    set(handles.ScaleEntry, 'String', '100');
end
scale=scale/100;

if scale==1
    bb=st{curfig}.bb;
    Diff=diff(abs(bb));
    
    Space=eye(4);
    if strcmpi(get(handles.LRButton, 'String'),'R')
        M=[ -1, 0, 0, -Diff(1);
            0, 1, 0, 0;
            0, 0, 1, 0;
            0, 0, 0, 1];
        Space=M*Space;
    end

    if strcmpi(get(handles.PAButton, 'String'),'P')
        M=[ 1, 0, 0, 0;
            0, -1, 0, -Diff(3);
            0, 0, 1, 0;
            0, 0, 0, 1];
        Space=M*Space;
    end
    if strcmpi(get(handles.SIButton, 'String'),'I')
        M=[ 1, 0, 0, 0;
            0, 1, 0, 0;
            0, 0, -1, -Diff(2);
            0, 0, 0, 1];
        Space=M*Space;
    end
    st{curfig}.Space=Space;
else
    Space=st{curfig}.Space;
    st{curfig}.Space(1:3,1:3)=Space(1:3,1:3)/(abs(Space(1:3,1:3))*scale);
    st{curfig}.Space(1:3,4)=st{curfig}.centre;
end
w_spm_orthviews('redraw');


% Hints: get(hObject,'String') returns contents of ScaleEntry as text
%        str2double(get(hObject,'String')) returns contents of ScaleEntry as a double

% --- Executes during object creation, after setting all properties.
function ScaleEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ScaleEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in CrosshairCheck.
function CrosshairCheck_Callback(hObject, eventdata, handles)
% hObject    handle to CrosshairCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
flag=get(handles.CrosshairCheck, 'Value');
if flag
    state='on';
else
    state='off';
end
w_spm_orthviews('xhairs', state)
% Hint: get(hObject,'Value') returns toggle state of CrosshairCheck

% --- Executes on button press in YokeCheck.
function YokeCheck_Callback(hObject, eventdata, handles)
% hObject    handle to YokeCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global st

curfig=handles.DPABI_fig;
curfig=w_Compatible2014bFig(curfig);

st{curfig}.yoke=1;
% Hint: get(hObject,'Value') returns toggle state of YokeCheck


% --- Executes on button press in NewButton.
function NewButton_Callback(hObject, eventdata, handles)
% hObject    handle to NewButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
New=SeeCAT_Viewer;
movegui(New, 'onscreen');

% --- Executes on button press in OverlayLabel.
function OverlayLabel_Callback(hObject, eventdata, handles)
% hObject    handle to OverlayLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DPABI_EGG
if isempty(DPABI_EGG)
    DPABI_EGG=0;
end
index=HeaderIndex(handles);
if ~index
    set(handles.OverlayLabel, 'Value', 0);
    DPABI_EGG=DPABI_EGG+1;
    if DPABI_EGG > 5
        msgbox('Don''t touch me, please!','modal');
        DPABI_EGG=0;
    end
    return
end

OverlayHeader=handles.OverlayHeaders{index};
RedrawOverlay(OverlayHeader, handles.DPABI_fig);
% Hint: get(hObject,'Value') returns toggle state of OverlayLabel


% --- Executes during object creation, after setting all properties.
function OverlayLabel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OverlayLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

function OverlayEntry_Callback(hObject, eventdata, handles)
% hObject    handle to OverlayEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global st;
% curfig=handles.DPABI_fig;
% curfig=w_Compatible2014bFig(curfig);
% 
% Num=get(handles.OverlayEntry, 'Value');
% if Num==1
%     set(handles.OverlayEntry, 'Value', st{curfig}.curblob+1);
%     return
% end
% set(handles.OverlayLabel, 'Value', 1);
% old_index=HeaderIndex(handles);
% OverlayHeader=handles.OverlayHeaders{old_index};
% RedrawOverlay(OverlayHeader, curfig);
% 
% st{curfig}.curblob=Num-1;
% index=HeaderIndex(handles, Num-1);
% OverlayHeader=handles.OverlayHeaders{index};
% RedrawOverlay(OverlayHeader, curfig);
% %state=w_OverlayConfig(handles.DPABI_fig, index, Num-1);
% 
% % handles=guidata(hObject);
% % if ~state
% %     set(handles.OverlayEntry, 'Value', st{curfig}.curblob+1);
% % end
% 
% guidata(hObject, handles)

% Hints: get(hObject,'String') returns contents of OverlayEntry as text
%        str2double(get(hObject,'String')) returns contents of OverlayEntry as a double


% --- Executes during object creation, after setting all properties.
function OverlayEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OverlayEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LeftButton.
function LeftButton_Callback(hObject, eventdata, handles)
% hObject    handle to LeftButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index=HeaderIndex(handles);
if ~index
    return
end
OverlayHeader=handles.OverlayHeaders{index};
oldTP=OverlayHeader.curTP;
curTP=oldTP-1;
set(handles.TimePointButton, 'String', sprintf('%d',curTP));
if curTP==1
    set(handles.LeftButton, 'Enable', 'Off');
end
set(handles.RightButton, 'Enable', 'On');

OverlayHeader=ChangeTP(OverlayHeader, curTP);
handles.OverlayHeaders{index}=OverlayHeader;
guidata(hObject, handles);
UpdateOverlayHandle(handles, 'On');

% --- Executes on button press in RightButton.
function RightButton_Callback(hObject, eventdata, handles)
% hObject    handle to RightButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index=HeaderIndex(handles);
if ~index
    return
end
OverlayHeader=handles.OverlayHeaders{index};
oldTP=OverlayHeader.curTP;
numTP=OverlayHeader.numTP;
curTP=oldTP+1;

set(handles.TimePointButton, 'String', sprintf('%d',curTP));
if curTP==numTP
    set(handles.RightButton, 'Enable', 'Off');
end
set(handles.LeftButton, 'Enable', 'On');

OverlayHeader=ChangeTP(OverlayHeader, curTP);
handles.OverlayHeaders{index}=OverlayHeader;
guidata(hObject, handles);
UpdateOverlayHandle(handles, 'On');

function TimePointButton_Callback(hObject, eventdata, handles)
% hObject    handle to TimePointButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index=HeaderIndex(handles);
if ~index
    return
end
OverlayHeader=handles.OverlayHeaders{index};
numTP=OverlayHeader.numTP;

curTP=str2double(get(handles.TimePointButton, 'String'));
if curTP<1 || curTP>numTP
    errordlg(sprintf('The number of time points: %d', numTP), 'Time Point Error');
    return;
end
if curTP==1
    set(handles.RightButton, 'Enable', 'On');
    set(handles.LeftButton, 'Enable', 'Off');
elseif curTP==numTP
    set(handles.RightButton, 'Enable', 'Off');
    set(handles.LeftButton, 'Enable', 'On');
else
    set(handles.RightButton, 'Enable', 'On');
    set(handles.LeftButton, 'Enable', 'On');    
end
OverlayHeader=ChangeTP(OverlayHeader, curTP);
handles.OverlayHeaders{index}=OverlayHeader;
guidata(hObject, handles);
UpdateOverlayHandle(handles, 'On');
% Hints: get(hObject,'String') returns contents of TimePointButton as text
%        str2double(get(hObject,'String')) returns contents of
%        TimePointButton as a double

% --- Executes during object creation, after setting all properties.
function TimePointButton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimePointButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function OverlayHeader=ChangeTP(OverlayHeader, curTP, curfig, curblob, OverlayVolumes)
global st

if nargin < 3
    curfig=gcf;
    curfig=w_Compatible2014bFig(curfig);
end
if nargin < 4
    curblob=st{curfig}.curblob;
end
if nargin < 5
    OverlayVolumes=[];
end

%cbarstring=OverlayHeader.cbarstring;
%NMax=OverlayHeader.NMax;
%NMin=OverlayHeader.NMin;
%PMin=OverlayHeader.PMin;
%PMax=OverlayHeader.PMax;
oldTP=OverlayHeader.curTP;
numTP=OverlayHeader.numTP;
if curTP~=oldTP
    if isempty(OverlayVolumes)
        OverlayFileName=OverlayHeader.fname;
        [OverlayVolume NewVox, NewHeader] = y_ReadRPI(OverlayFileName, curTP);
    else
        OverlayVolume=OverlayVolumes(:,:,:,curTP);
    end
        
    OverlayHeader.Raw = OverlayVolume;
else
    return;
end
NMax=min(OverlayVolume(:));
PMax=max(OverlayVolume(:));
OverlayHeader.NMax=NMax;
OverlayHeader.PMax=PMax;
OverlayHeader.curTP=curTP;
OverlayHeader.numTP=numTP;
OverlayHeader.IsSelected=1;

OverlayHeader=RedrawOverlay(OverlayHeader, curfig);

% --- Executes on button press in OverlayButton.
function OverlayButton_Callback(hObject, eventdata, handles)
% hObject    handle to OverlayButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[File , Path]=uigetfile({'*.img;*.nii;*.nii.gz','Brain Image Files (*.img;*.nii;*.nii.gz)';'*.*', 'All Files (*.*)';}, ...
    'Pick Overlay File' , pwd); 

if ~ischar(File)
    return;
end
OverlayFile=fullfile(Path, File);
set(handles.OverlayEntry, 'String', File, 'Enable', 'On');
handles.OverlayFileName=cell(1, 1);
handles.OverlayFileName{1}=OverlayFile;
%Added by YAN Chao-Gan, 131005. To handle .nii.gz files.
FileNameTemp = OverlayFile;
[pathstr, name, ext] = fileparts(FileNameTemp);
if isempty(ext)
    FileNameTemp = fullfile(pathstr,[name '.nii']);
    if ~exist(FileNameTemp,'file')
        FileNameTemp = fullfile(pathstr,[name '.nii.gz']);
        [pathstr, name, ext] = fileparts(FileNameTemp);
    end
    if ~exist(FileNameTemp,'file')
        FileNameTemp = fullfile(pathstr,[name '.hdr']);
    end
end
if ~exist(FileNameTemp,'file')
    error(['File doesn''t exist: ',fullfile(pathstr,[name ext])]);
end
if strcmpi(ext,'.gz')
    gunzip(FileNameTemp);
    FileNameTemp = fullfile(pathstr,[name]);
end

Nii = nifti(FileNameTemp);
numTP = size(Nii.dat,4);
curTP = 1;
if numTP>1
    State='On';
    stringTP=num2str(curTP);
else
    State='Off';
    stringTP='';
end
set(handles.LeftButton, 'Enable', 'Off');
set(handles.RightButton, 'Enable', State);
set(handles.TimePointButton, 'Enable', State);
set(handles.TimePointButton, 'String', stringTP);

if strcmpi(ext,'.gz')
    delete(FileNameTemp);
end
%END%Added by YAN Chao-Gan, 131005. To handle .nii.gz files.

[OverlayVolume OverlayVox OverlayHeader] = y_ReadRPI(OverlayFile, curTP);
OverlayHeader=w_ExtendHeader(OverlayHeader);
OverlayHeader.Vox=OverlayVox;

PN_Flag='';
NMax=min(OverlayVolume(:));
if NMax >= 0
    PN_Flag='+';
    NMax=0;
end
NMin = 0;
PMin = 0;
PMax=max(OverlayVolume(:));
if PMax <= 0;
    PN_Flag='-';
    PMax=0;
end

OverlayHeader.IsSelected=1;
OverlayHeader.Raw=OverlayVolume;
OverlayHeader.Data=OverlayVolume;
OverlayHeader.NMax=NMax;
OverlayHeader.NMin=NMin;
OverlayHeader.PMin=PMin;
OverlayHeader.PMax=PMax;
OverlayHeader.cbarstring='0';
OverlayHeader.numTP=numTP;
OverlayHeader.curTP=curTP;
OverlayHeader.PN_Flag=PN_Flag;
OverlayHeader.Alpha=1;
OverlayHeader.ColorMap=colormap('jet(64)');
OverlayHeader.ColorMapIndex=1;

SendHeader=OverlayHeader;
if get(handles.OnlyPos, 'Value')
    SendHeader.Data(SendHeader.Data < 0)=0;
end
if get(handles.OnlyNeg, 'Value')
    SendHeader.Data(SendHeader.Data > 0)=0;
end
if get(handles.OnlyUnder, 'Value')
    SendHeader.Data(SendHeader.Data~= 0)=0;
end
SendHeader=SetCSize(SendHeader);

ColorMap = AdjustColorMap(OverlayHeader.ColorMap,...
    [0.75 0.75 0.75],...
    NMax,...
    NMin,...
    PMin,...
    PMax,...
    OverlayHeader.PN_Flag);
w_spm_orthviews('removeblobs', handles.DPABI_fig);
w_spm_orthviews('Addtruecolourimage',...
    w_Compatible2014bFig(handles.DPABI_fig),...
    SendHeader,...
    ColorMap,...
    SendHeader.Alpha,...
    PMax,...
    NMax);
w_spm_orthviews('redrawcolourbar', w_Compatible2014bFig(handles.DPABI_fig), 1);
w_spm_orthviews('Redraw', w_Compatible2014bFig(handles.DPABI_fig));

handles.OverlayHeaders=cell(5, 1);
handles.OverlayHeaders{1}=OverlayHeader;
guidata(hObject, handles);
UpdateOverlayHandle(handles, 'On');

function UpdateOverlayHandle(handles, State)
index=HeaderIndex(handles);
OverlayHeader=handles.OverlayHeaders{index};
%% Overlay Color
ColorMapIndex=OverlayHeader.ColorMapIndex;
set(handles.ColormapPopup, 'Value', ColorMapIndex, 'Enable', State);
PN_Flag=OverlayHeader.PN_Flag;
if isempty(PN_Flag) % Pos/Neg
    ColorbarType=3;
elseif strcmpi(PN_Flag, '+') % Full Pos
    ColorbarType=1;
else
    ColorbarType=2;
end
set(handles.ColorbarTypePopup, 'Value', ColorbarType, 'Enable', State);
Alpha=OverlayHeader.Alpha;
set(handles.AlphaSlider, 'Value', Alpha, 'Enable', State);
set(handles.AlphaEntry, 'String', num2str(Alpha), 'Enable', State);
%% Overlay Threshold
PMin=OverlayHeader.PMin;
NMin=OverlayHeader.NMin;
if abs(PMin)>abs(NMin)
    Thrd=PMin;
else
    Thrd=abs(NMin);
end
PMax=OverlayHeader.PMax;
NMax=OverlayHeader.NMax;
if abs(PMax)>abs(NMax)
    Max=PMax;
else
    Max=abs(NMax);
end
set(handles.ThrdEntry, 'String', num2str(Thrd));
if Thrd/Max>1
    warning('Threshold Too High!');
    set(handles.ThrdSlider, 'Value', 1);
else
    set(handles.ThrdSlider, 'Value', Thrd/Max);
end
%% Overlay CS
RMM=OverlayHeader.RMM;
switch RMM
    case 6
        CSType=1;
    case 18
        CSType=2;
    case 26
        CSType=3;
end
set(handles.CSTypePopup, 'Value', CSType, 'Enable', State);
CSize=OverlayHeader.CSize;
V=prod(OverlayHeader.Vox);
CVSize=CSize/V;
set(handles.CSEntry, 'String', CVSize, 'Enable', State);

%Correct
set(handles.CorrectPopup, 'Enable', State);

%Thres
set(handles.OnlyPos, 'Enable', State);
set(handles.OnlyNeg, 'Enable', State);
set(handles.OnlyUnder, 'Enable', State);
set(handles.ThrdEntry, 'Enable', State);
set(handles.ThrdSlider, 'Enable', State);
set(handles.PEntry, 'Enable', State);
set(handles.DfButton, 'Enable', State);
%Report
set(handles.FindPeakBtn, 'Enable', State);
set(handles.CIReportBtn, 'Enable', State);

function UnderlayEntry_Callback(hObject, eventdata, handles)
% hObject    handle to UnderlayEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.UnderlayFileName)
    [DPABIPath, fileN, extn] = fileparts(which('SeeCAT.m'));
    TemplatePath=fullfile(DPABIPath, 'templates');
    
    UnderlayFileName=fullfile(TemplatePath,'ch2.nii');
    [Path, Name, Ext]=fileparts(UnderlayFileName);
else
    [Path, Name, Ext]=fileparts(handles.UnderlayFileName);
end
NewName=get(handles.UnderlayEntry, 'String');
[NewPath, NewName, NewExt]=fileparts(NewName);
if isempty(NewPath)
    NewPath=Path;
end
NewFileName=[NewPath, filesep, NewName, NewExt];
if ~exist(NewFileName, 'file')
    errordlg('Image File not found', 'File Error');
    set(handles.UnderlayEntry, 'String', [Name, Ext]);
else
    handles.UnderlayFileName=NewFileName;
    guidata(hObject, handles);    
    ShowUnderlay(handles);
end
% Hints: get(hObject,'String') returns contents of UnderlayEntry as text
%        str2double(get(hObject,'String')) returns contents of UnderlayEntry as a double


% --- Executes during object creation, after setting all properties.
function UnderlayEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UnderlayEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in UnderlayButton.
function UnderlayButton_Callback(hObject, eventdata, handles)
% hObject    handle to UnderlayButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.UnderlayFileName)
    [File , Path]=uigetfile({'*.img;*.nii;*.nii.gz','Brain Image Files (*.img;*.nii;*.nii.gz)';'*.*', 'All Files (*.*)';}, ...
        'Pick Underlay File' , pwd);
else    
    [Path, Name, Ext]=fileparts(handles.UnderlayFileName);
    [File , Path]=uigetfile({'*.img;*.nii;*.nii.gz','Brain Image Files (*.img;*.nii;*.nii.gz)';'*.*', 'All Files (*.*)';}, ...
        'Pick Underlay File' , [Path , filesep , Name , Ext]);
end
if ~ischar(File)
    return;
end
UnderlayFileName=[Path, File];
handles.UnderlayFileName=UnderlayFileName;
handles.UserDefinedFileName=UnderlayFileName;
set(handles.UnderlayEntry, 'String', File);
guidata(hObject, handles);
ShowUnderlay(handles);

function ShowUnderlay(handles)
if isempty(handles.UnderlayFileName)
    [DPABIPath, fileN, extn] = fileparts(which('SeeCAT.m'));
    TemplatePath=fullfile(DPABIPath, 'templates');
    UnderlayFileName=fullfile(TemplatePath,'ch2.nii');
else
    UnderlayFileName=handles.UnderlayFileName;
end
[UnderlayVolume UnderlayVox UnderlayHeader] = y_ReadRPI(UnderlayFileName);
UnderlayHeader.Data = UnderlayVolume;
UnderlayHeader.Vox  = UnderlayVox;

w_spm_orthviews('Image',UnderlayHeader);

Max=length(get(handles.TemplatePopup, 'String'));
set(handles.TemplatePopup, 'Value', Max-1);

function handles=NoneUnderlay(handles)
global st
curfig=handles.DPABI_fig;
curfig=w_Compatible2014bFig(curfig);
index=HeaderIndex(handles);
if ~index
    st{curfig}.vols{1}.Data=zeros(st{curfig}.vols{1}.dim);
    w_spm_orthviews('Redraw');
else
    OverlayHeader=handles.OverlayHeaders{index};
    OverlayHeader.Data=zeros(OverlayHeader.dim);
    w_spm_orthviews('Image', OverlayHeader);
    handles.UnderlayFileName=OverlayHeader.fname;
end


% --- Executes on selection change in TemplatePopup.
function TemplatePopup_Callback(hObject, eventdata, handles)
% hObject    handle to TemplatePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
flag=get(handles.TemplatePopup, 'Value');
Max=length(get(handles.TemplatePopup, 'String'));
[DPABIPath, fileN, extn] = fileparts(which('SeeCAT.m'));
switch flag
    case 1
        File='ch2.nii';
    case 2
        File='ch2bet.nii';
    case 3
        File='ch2better.nii';
    case 4
        File='mni_icbm152_t1_tal_nlin_asym_09c.nii';
    case 5
        UnderlayFileName=handles.UserDefinedFileName;
        [UnderlayVolume UnderlayVox UnderlayHeader] = y_ReadRPI(UnderlayFileName);
        UnderlayHeader.Data = UnderlayVolume;
        UnderlayHeader.Vox  = UnderlayVox;
        w_spm_orthviews('Image',UnderlayHeader);
        return;
    case Max
        handles=NoneUnderlay(handles);
        set(handles.UnderlayEntry, 'String', '');
        guidata(hObject, handles);
        return;
    otherwise
        return;
end
UnderlayFileName=[DPABIPath,filesep,'templates',filesep, File];
[UnderlayVolume UnderlayVox UnderlayHeader] = y_ReadRPI(UnderlayFileName);
UnderlayHeader.Data = UnderlayVolume;
UnderlayHeader.Vox  = UnderlayVox;
w_spm_orthviews('Image',UnderlayHeader);

set(handles.UnderlayEntry, 'String', File);
set(handles.TemplatePopup, 'Value', flag);

handles.UnderlayFileName=UnderlayFileName;
guidata(hObject, handles);

% Hints: contents = get(hObject,'String') returns TemplatePopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TemplatePopup


% --- Executes during object creation, after setting all properties.
function TemplatePopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TemplatePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function ThrdSlider_Callback(hObject, eventdata, handles)
% hObject    handle to ThrdSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ThrdListener_Callback(hObject, eventdata, hObject);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

function ThrdListener_Callback(hObject, eventdata, hFig)
handles=guidata(hFig);
curfig=handles.DPABI_fig;
curfig=w_Compatible2014bFig(curfig);

value=get(handles.ThrdSlider, 'Value');
index=HeaderIndex(handles);
if ~index
    return
end

OverlayHeader=handles.OverlayHeaders{index};

PMax=OverlayHeader.PMax;
NMax=OverlayHeader.NMax;
Max=PMax;

if abs(NMax) > abs(PMax)
    Max=abs(NMax);
end

Thrd=Max*value;

set(handles.ThrdEntry, 'String', sprintf('%g', Thrd));
OverlayHeader=ChangeThrd(OverlayHeader, -Thrd, Thrd, curfig);
handles.OverlayHeaders{index}=OverlayHeader;
guidata(curfig, handles);
UpdateOverlayHandle(handles, 'On');

switch upper(OverlayHeader.TestFlag)
    case 'T'
        Df=OverlayHeader.Df;
        P=2*(1-tcdf(abs(Thrd), Df));
    case 'R'
        Df=OverlayHeader.Df;
        T=sqrt(Df*(Thrd^2/(1-Thrd^2)));
        P=2*(1-tcdf(abs(Thrd)*T, Df));
    case 'F'
        Df1=OverlayHeader.Df;
        Df2=OverlayHeader.Df2;
        P=(1-fcdf(Thrd, Df1, Df2));
    case 'Z'
        P=2*(1-normcdf(abs(Thrd)));
    otherwise
        return;
end
set(handles.PEntry, 'String', num2str(P));

function [OverlayHeader, SendHeader]=RedrawOverlay(SendHeader, curfig, curblob)
global st
if nargin<2
    curfig=gcf;
end
curfig=w_Compatible2014bFig(curfig);

MainHandle=guidata(curfig);
if nargin<3
    curblob=st{curfig}.curblob;
end
Alpha=SendHeader.Alpha;
PMax=SendHeader.PMax;
PMin=SendHeader.PMin;
NMin=SendHeader.NMin;
NMax=SendHeader.NMax;
PN_Flag=SendHeader.PN_Flag;
ColorMap=SendHeader.ColorMap;

OverlayVolume=SendHeader.Raw;
OverlayVolume(~SendHeader.Mask)=0;
if NMax >= 0
    OverlayVolume(OverlayVolume<0) = 0;
end
if PMax <= 0
    OverlayVolume(OverlayVolume>0) = 0;
end

SendHeader.Data=OverlayVolume.*...
    ((OverlayVolume < NMin) + (OverlayVolume > PMin));
OverlayHeader=SendHeader;

if SendHeader.CSize
    SendHeader=SetCSize(SendHeader);
end
if ~isempty(SendHeader.ThrdIndex)
    Tokens=regexp(SendHeader.ThrdIndex, '(\d+\.*\d*:?\d*\.*\d*)', 'tokens');
    if ~isempty(Tokens)
        L=false(size(SendHeader.Data));
        for t=1:numel(Tokens)
            Num=Tokens{t};
            Num=str2num(Num{1});
            LL=false(size(SendHeader.Data));
            if length(Num)==1
                LL(SendHeader.Data==Num)=1;
            else
                LL(SendHeader.Data>=Num(1) & SendHeader.Data<=Num(end))=1;
            end
            L=L+LL;
        end
        SendHeader.Data(~L)=0;
    end
end
%Only +/-/Display Current Overlay/No Overlay
if get(MainHandle.OnlyPos, 'Value')
    SendHeader.Data(OverlayVolume < 0) = 0;
elseif get(MainHandle.OnlyNeg, 'Value')
    SendHeader.Data(OverlayVolume > 0) = 0;
elseif get(MainHandle.OnlyUnder, 'Value')
    SendHeader.Data(OverlayVolume~= 0) = 0;
end

ColorMap = AdjustColorMap(ColorMap,...
    [0.75 0.75 0.75],...
    NMax,...
    NMin,...
    PMin,...
    PMax,...
    PN_Flag);

w_spm_orthviews('Settruecolourimage',...
    curfig,...
    SendHeader,...
    ColorMap,...
    Alpha,...
    PMax,...
    NMax,...
    curblob);

handles=guidata(curfig);
if get(handles.TemplatePopup, 'Value')==length(get(handles.TemplatePopup, 'String'))
    handles=NoneUnderlay(handles);
    guidata(curfig, handles);
else
    w_spm_orthviews('Redraw', curfig);
end

function OverlayHeader=ChangeThrd(OverlayHeader, NMin, PMin, curfig)
if nargin<4
    curfig=gcf;
    curfig=w_Compatible2014bFig(curfig);
end
OverlayHeader.NMin=NMin;
OverlayHeader.PMin=PMin;
OverlayHeader=RedrawOverlay(OverlayHeader, curfig);

% --- Executes during object creation, after setting all properties.
function ThrdSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ThrdSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function ThrdEntry_Callback(hObject, eventdata, handles)
% hObject    handle to ThrdEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index=HeaderIndex(handles);
if ~index
    set(handles.ThrdEntry, 'String', []);
    return
end

Thrd=str2double(get(handles.ThrdEntry, 'String'));
if isnan(Thrd)
    set(handles.ThrdEntry, 'String', []);
    return;
end

curfig=handles.DPABI_fig;
curfig=w_Compatible2014bFig(curfig);

OverlayHeader=handles.OverlayHeaders{index};

PMax=OverlayHeader.PMax;
NMax=OverlayHeader.NMax;
Max=PMax;

if abs(NMax) > abs(PMax)
    Max=abs(NMax);
end
Value=Thrd/Max;

set(handles.ThrdSlider, 'Value', Value);
set(handles.ThrdEntry, 'String', sprintf('%g', Thrd));
OverlayHeader=ChangeThrd(OverlayHeader, -Thrd, Thrd, curfig);
handles.OverlayHeaders{index}=OverlayHeader;
guidata(curfig, handles);
% Hints: get(hObject,'String') returns contents of ThrdEntry as text
%        str2double(get(hObject,'String')) returns contents of ThrdEntry as a double


% --- Executes during object creation, after setting all properties.
function ThrdEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ThrdEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function PEntry_Callback(hObject, eventdata, handles)
% hObject    handle to PEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index=HeaderIndex(handles);
if ~index
    set(handles.PEntry, 'String', []);
    return
end
OverlayHeader=handles.OverlayHeaders{index};
P=str2double(get(handles.PEntry, 'String'));
if isnan(P)
    set(handles.PEntry, 'String', []);
    return;
end
switch upper(OverlayHeader.TestFlag)
    case 'T'
        Df=OverlayHeader.Df;
        Thrd=tinv(1-P/2, Df);
    case 'R'
        Df=OverlayHeader.Df;
        T=tinv(1-P/2, Df);
        Thrd=sqrt(T^2/(Df+T^2));
    case 'F'
        Df1=OverlayHeader.Df;
        Df2=OverlayHeader.Df2;
        Thrd=finv(1-P, Df1, Df2);
    case 'Z'
        Thrd=norminv(1-P/2);
    otherwise
        return;
end
set(handles.ThrdEntry, 'String', num2str(Thrd));
%Max=max([abs(OverlayHeader.PMax),abs(OverlayHeader)]);
%Max
OverlayHeader=ChangeThrd(OverlayHeader, -Thrd, Thrd, handles.DPABI_fig);
handles.OverlayHeaders{index}=OverlayHeader;
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of PEntry as text
%        str2double(get(hObject,'String')) returns contents of PEntry as a double


% --- Executes during object creation, after setting all properties.
function PEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in DfButton.
function DfButton_Callback(hObject, eventdata, handles)
% hObject    handle to DfButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index=HeaderIndex(handles);
if ~index
    return
end
OverlayHeader=handles.OverlayHeaders{index};
Para=w_ChangeDf(OverlayHeader);
if ~iscell(Para)
    return;
end
OverlayHeader.TestFlag=Para{1};
OverlayHeader.Df=Para{2};
OverlayHeader.Df2=Para{3};

handles.OverlayHeaders{index}=OverlayHeader;
guidata(hObject, handles);

function PosNegNon(OverlayHeaders, flag, curfig)
global st
if nargin<3
    curfig=gcf;
    curfig=w_Compatible2014bFig(curfig);
end

if isfield(st{curfig}.vols{1}, 'blobs')
    blob=0;
    for i=1:numel(OverlayHeaders)
        if ~isempty(OverlayHeaders{i}) && ...
                OverlayHeaders{i}.IsSelected
            blob=blob+1;
            OverlayHeader=OverlayHeaders{i};
            if strcmpi(flag, '+')
                OverlayHeader.Data(OverlayHeader.Data < 0)=0;
            elseif strcmpi(flag, '-')
                OverlayHeader.Data(OverlayHeader.Data > 0)=0;
            elseif strcmpi(flag, 'N')
                OverlayHeader.Data(OverlayHeader.Data~= 0)=0;
            end
            OverlayHeader=SetMask(OverlayHeader);
            OverlayHeader=SetCSize(OverlayHeader);
            st{curfig}.vols{1}.blobs{blob}.vol=OverlayHeader;
        else
            continue;
        end
    end
end

% --- Executes on button press in OnlyPos.
function OnlyPos_Callback(hObject, eventdata, handles)
% hObject    handle to OnlyPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
curfig=handles.DPABI_fig;
curfig=w_Compatible2014bFig(curfig);

flag=get(handles.OnlyPos, 'Value');

if flag
    set(handles.OnlyNeg,   'Value', ~flag);
    set(handles.OnlyUnder, 'Value', ~flag);
    PosNegNon(handles.OverlayHeaders, '+', curfig);
else
    PosNegNon(handles.OverlayHeaders, 'A', curfig);
end

w_spm_orthviews('Redraw');
guidata(hObject, handles);

% Hint: get(hObject,'Value') returns toggle state of OnlyPos


% --- Executes on button press in OnlyNeg.
function OnlyNeg_Callback(hObject, eventdata, handles)
% hObject    handle to OnlyNeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
curfig=handles.DPABI_fig;
curfig=w_Compatible2014bFig(curfig);

flag=get(handles.OnlyNeg, 'Value');

if flag
    set(handles.OnlyPos,   'Value', ~flag);
    set(handles.OnlyUnder, 'Value', ~flag);
    PosNegNon(handles.OverlayHeaders, '-', curfig);
else
    PosNegNon(handles.OverlayHeaders, 'A', curfig);
end

w_spm_orthviews('Redraw');
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of OnlyNeg


% --- Executes on button press in OnlyUnder.
function OnlyUnder_Callback(hObject, eventdata, handles)
% hObject    handle to OnlyUnder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
curfig=handles.DPABI_fig;
curfig=w_Compatible2014bFig(curfig);

flag=get(handles.OnlyUnder, 'Value');

if flag
    set(handles.OnlyPos,   'Value', ~flag);
    set(handles.OnlyNeg, 'Value', ~flag);
    PosNegNon(handles.OverlayHeaders, 'N', curfig);
else
    PosNegNon(handles.OverlayHeaders, 'A', curfig);
end

w_spm_orthviews('Redraw');
guidata(hObject, handles);

% Hint: get(hObject,'Value') returns toggle state of OnlyUnder


% --- Executes on selection change in ClusterPopup.
function ClusterPopup_Callback(hObject, eventdata, handles)
% hObject    handle to ClusterPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index=HeaderIndex(handles);
if ~index
    return
end
Value=get(handles.ClusterPopup, 'Value');
OverlayHeader=handles.OverlayHeaders{index};
switch Value
    case 1
        return
    case 2 %Set Cluster Size
        OverlayHeader=w_ClusterSize(OverlayHeader);
        if isempty(OverlayHeader)
            return
        end
        OverlayHeader=RedrawOverlay(OverlayHeader);
        handles.OverlayHeaders{index}=OverlayHeader;
    case 3 %Save Single Cluster
        [File , Path]=uiputfile({'*.img;*.nii;*.nii.gz','Brain Image Files (*.img;*.nii;*.nii.gz)';'*.*', 'All Files (*.*)';}, ...
            'Save All Clusters' , OverlayHeader.fname);
        if ~ischar(File)
            return
        end
        [OverlayHeader, SendHeader]=RedrawOverlay(OverlayHeader);
        
        [Path, FileName, Ext]=fileparts(fullfile(Path, File));
        MaskName=sprintf('%s_mask',FileName);
        Data=SingleCluster(SendHeader);
        if isempty(Data)
            return
        end
        Mask=logical(Data);
        y_Write(Data, SendHeader, fullfile(Path, [FileName, Ext]));
        y_Write(Mask, SendHeader, fullfile(Path, [MaskName, Ext]));
    case 4 %Save All Clusters
        [File , Path]=uiputfile({'*.img;*.nii;*.nii.gz','Brain Image Files (*.img;*.nii;*.nii.gz)';'*.*', 'All Files (*.*)';}, ...
            'Save All Clusters' , OverlayHeader.fname);
        if ~ischar(File)
            return
        end
        [OverlayHeader, SendHeader]=RedrawOverlay(OverlayHeader);
        
        [Path, FileName, Ext]=fileparts(fullfile(Path, File));
        MaskName=sprintf('%s_mask',FileName);
        Data=SendHeader.Data;
        Mask=logical(Data);
        y_Write(Data, SendHeader, fullfile(Path, [FileName, Ext]));
        y_Write(Mask, SendHeader, fullfile(Path, [MaskName, Ext]));
    case 5 %Find Peak in this Cluster
        OverlayHeader=SetCSize(OverlayHeader);
        OverlayHeader=SetPosNeg(OverlayHeader, handles);
        FindPeak(OverlayHeader);
    case 6 %Volume Percentage
        OverlayHeader=w_Percentage(OverlayHeader);
        if isempty(OverlayHeader)
            return
        end
        OverlayHeader=RedrawOverlay(OverlayHeader);
        handles.OverlayHeaders{index}=OverlayHeader;        
    case 7 %FDR
        OverlayHeader=w_FDRCorrection(OverlayHeader);
        if isempty(OverlayHeader)
            return
        end
        OverlayHeader=RedrawOverlay(OverlayHeader, handles.DPABI_fig);
        handles.OverlayHeaders{index}=OverlayHeader;        
    case 8 %GRF
        OverlayHeader=w_GRFCorrection(OverlayHeader);
        if isempty(OverlayHeader)
            return
        end
        [OverlayHeader, SendHeader]=RedrawOverlay(OverlayHeader);
        %OverlayHeader.Data = SendHeader.Data; %OverlayHeader=SetCSize(OverlayHeader); %YAN Chao-Gan, 140822. Need to save the data after setting cluster size.
        handles.OverlayHeaders{index}=OverlayHeader;        
    case 9 %AlphaSim
        w_AlphaSimCorrection(OverlayHeader);
    case 10 %Cluster Report
        OverlayHeader=SetCSize(OverlayHeader);               
        OverlayHeader=SetPosNeg(OverlayHeader, handles); %Add by Sandy to set cluster size at first when someone use Cluster Report       
        y_ClusterReport(OverlayHeader.Data, OverlayHeader, OverlayHeader.RMM);
end
guidata(hObject, handles);
% Hints: contents = get(hObject,'String') returns ClusterPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ClusterPopup


% --- Executes during object creation, after setting all properties.
function ClusterPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ClusterPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in MorePopup.
function MorePopup_Callback(hObject, eventdata, handles)
% hObject    handle to MorePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global st
curfig=handles.DPABI_fig;
curfig=w_Compatible2014bFig(curfig);
curblob=st{curfig}.curblob;

Value=get(handles.MorePopup, 'Value');
switch Value
    case 1
        return
    case 2 %Time Course
        Headers=cell(numel(handles.OverlayHeaders), 1);
        f=0;
        for i=1:numel(handles.OverlayHeaders)
            if ~isempty(handles.OverlayHeaders{i}) && ...
                    handles.OverlayHeaders{i}.IsSelected && ...
                    handles.OverlayHeaders{i}.numTP > 1
                Headers{i}=handles.OverlayHeaders{i};
                f=1;
            end
        end
        if ~f
            return;
        end
        w_TimeCourse(handles.DPABI_fig, Headers);
    case 3 %Call BrainNet Viewer
        if ~(exist('BrainNet.m'))
            msgbox('The surface view is based on Mingrui Xia''s BrainNet Viewer. Please install BrainNet Viewer 1.1 or later version at first (http://www.nitrc.org/projects/bnv/).','REST Slice Viewer', 'modal');
            return
        end
        
        index=HeaderIndex(handles);
        OverlayHeader=handles.OverlayHeaders{index};
        OverlayHeader=SetPosNeg(OverlayHeader, handles);
        OverlayHeader=SetCSize(OverlayHeader);
        CVSize=OverlayHeader.CSize/prod(OverlayHeader.Vox);
        
        cbarstring = OverlayHeader.cbarstring;
        if cbarstring(end)=='+' || cbarstring(end)=='-'
            PN_Flag=cbarstring(end);
            cbarstring=cbarstring(1:end-1);
        else
            PN_Flag=[];
        end
        cbar=str2double(cbarstring);
        
        if cbar==0 && isfield(SendHeader, 'ColorMap')
            ColorMap=OverlayHeader.ColorMap;
        else
            if isnan(cbar)
                ColorMap = colormap(cbarstring);
            else
                ColorMap = y_AFNI_ColorMap(cbar);
            end
        end

%         ColorMap = AdjustColorMap(ColorMap,...
%             [0.75 0.75 0.75],...
%             OverlayHeader.NMax,...
%             OverlayHeader.NMin,...
%             OverlayHeader.PMin,...
%             OverlayHeader.PMax);

        [BrainNetViewerPath, fileN, extn] = fileparts(which('BrainNet.m'));
        SurfFileName=[BrainNetViewerPath,filesep,'Data',filesep,'SurfTemplate',filesep,'BrainMesh_ICBM152_smoothed.nv'];
        y_CallBrainNetViewer(OverlayHeader.Data,...
            OverlayHeader.NMin, OverlayHeader.PMin,...
            CVSize, OverlayHeader.RMM,...
            SurfFileName, 'MediumView', ColorMap,...
            OverlayHeader.NMax, OverlayHeader.PMax,...
            OverlayHeader,...
            PN_Flag);
    case 4 % Set Range
        index=HeaderIndex(handles);
        if ~index
            return
        end
        OverlayHeader=handles.OverlayHeaders{index};
        ThrdIndex=OverlayHeader.ThrdIndex;
        ThrdIndex=inputdlg('Set threshold range (e.g., "2:5,10,12" means all the values within 2 and 5 (including decimals) and values of 10 and 12.',...
            'Range of Overlay''s Thrdhold',...
            1, {ThrdIndex});
        if isempty(ThrdIndex)
            return
        end
        OverlayHeader.ThrdIndex=ThrdIndex{1};
        OverlayHeader=RedrawOverlay(OverlayHeader);
        handles.OverlayHeaders{index}=OverlayHeader;
    case 5 % Set Colorbar Manually
        index=HeaderIndex(handles);
        if ~index
            return
        end
        OverlayHeader=handles.OverlayHeaders{index};
        
        colormap(OverlayHeader.ColorMap);
        w_colormapeditor(handles.DPABI_fig);
    case 6 % Transparency
        Transparency=1-st{curfig}.vols{1}.blobs{curblob}.colour.prop;
        Transparency=inputdlg('Set Transparency (Default = 0.2)',...
            'Transparency of Overlay',...
            1, {num2str(Transparency)});
        if isempty(Transparency)
            return
        end
        st{curfig}.vols{1}.blobs{curblob}.colour.prop=1-str2double(Transparency{1});
        w_spm_orthviews('Redraw', curfig);
    case 7 %Save as Picture
        [File, Path] = uiputfile({'*.tiff';'*.jpeg';'*.png';'*.bmp'},...
            'Save Picture As');
        if ~ischar(File)
            return;
        end
        [Path, Name, Ext]=fileparts(fullfile(Path, File));
        TName=sprintf('%s_Transverse', Name);
        CName=sprintf('%s_Coronal', Name);
        SName=sprintf('%s_Sagittal', Name);
        TData=getframe(st{curfig}.vols{1}.ax{1}.ax);
        CData=getframe(st{curfig}.vols{1}.ax{2}.ax);
        SData=getframe(st{curfig}.vols{1}.ax{3}.ax);
 
        saveas(curfig, fullfile(Path, [Name, Ext]));
        eval(['print -r300 -dtiff -noui ''',fullfile(Path, [Name, '_300dpi', Ext]),''';']); %YAN Chao-Gan, 140806.
        
        imwrite(TData.cdata, fullfile(Path, [TName, Ext]));
        imwrite(CData.cdata, fullfile(Path, [CName, Ext]));
        imwrite(SData.cdata, fullfile(Path, [SName, Ext]));
    case 8 %Save colorbar as Picture
        [File, Path] = uiputfile({'*.tiff';'*.jpeg';'*.png';'*.bmp'},...
            'Save Picture As');
        if ~ischar(File)
            return;
        end
        
        index=HeaderIndex(handles);
        if ~index
            return
        end        
        OverlayHeader=handles.OverlayHeaders{index};
        
        cbarstring=OverlayHeader.cbarstring;
        if cbarstring(end)=='+' || cbarstring(end)=='-'
            PN_Flag=cbarstring(end);
            cbarstring=cbarstring(1:end-1);
        else
            PN_Flag=[];
        end        
        cbar=str2double(cbarstring);
        if cbar==0 && isfield(OverlayHeader, 'ColorMap')
            ColorMap=OverlayHeader.ColorMap;
        else
            if isnan(cbar)
                ColorMap = colormap(cbarstring);
            else
                ColorMap = y_AFNI_ColorMap(cbar);
            end
        end
        
        %ColorMap=flipdim(ColorMap, 1);
        L=size(ColorMap, 1);
        NMax=OverlayHeader.NMax;
        NMin=OverlayHeader.NMin;
        PMin=OverlayHeader.PMin;
        PMax=OverlayHeader.PMax;
        if NMax==0 && NMax==NMin;
            Scale=(L:-1:L/2+1)';
        elseif PMax==0 && PMax==PMin
            Scale=(L/2:-1:1)';
        else
            Scale=(1:L)';
        end
        
        imwrite(imresize(Scale, [320*L, 200]),...
            ColorMap, fullfile(Path, File));
    case 9 %Save colorbar as MAT
        index=HeaderIndex(handles);
        if ~index
            return
        end
        
        [File, Path] = uiputfile('Colormap.mat',...
            'Save Colormap As MAT');
        if ~ischar(File)
            return;
        end
                
        OverlayHeader=handles.OverlayHeaders{index};
        
        cbarstring=OverlayHeader.cbarstring;
        if cbarstring(end)=='+' || cbarstring(end)=='-'
            PN_Flag=cbarstring(end);
            cbarstring=cbarstring(1:end-1);
        else
            PN_Flag=[];
        end        
        cbar=str2double(cbarstring);
        
        if cbar==0 && isfield(OverlayHeader, 'ColorMap')
            ColorMap=OverlayHeader.ColorMap;
        else
            if isnan(cbar)
                ColorMap = colormap(cbarstring);
            else
                ColorMap = y_AFNI_ColorMap(cbar);
            end
        end
        OverlayHeader.ColorMap=ColorMap;
        save(fullfile(Path, File), 'ColorMap');
    case 10 %Load colorbar from MAT
        index=HeaderIndex(handles);
        if ~index
            return
        end     
        
        [File, Path] = uigetfile('*.mat',...
            'Load Colormap from MAT');
        if ~ischar(File)
            return;
        end
         
        Temp=load(fullfile(Path, File));
        ColorMap=Temp.ColorMap;
        OverlayHeader=handles.OverlayHeaders{index};
        OverlayHeader.ColorMap=ColorMap;
        OverlayHeader.cbarstring='0';
        OverlayHeader=RedrawOverlay(OverlayHeader, handles.DPABI_fig);
        
        handles.OverlayHeader{index}=OverlayHeader;
end
guidata(hObject, handles);
% Hints: contents = get(hObject,'String') returns MorePopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MorePopup


% --- Executes during object creation, after setting all properties.
function MorePopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MorePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in SIButton.
function SIButton_Callback(hObject, eventdata, handles)
% hObject    handle to SIButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global st
curfig=handles.DPABI_fig;
curfig=w_Compatible2014bFig(curfig);
SI={'S','I'};
state=get(handles.SIButton, 'String');
state=strcmpi(state, 'S')+1;
state=SI{state};

bb=st{curfig}.bb;
offset=-bb(1,2)-bb(2,2);
M=[ 1, 0, 0, 0;
    0, 1, 0, 0;
    0, 0, -1, offset;
    0, 0, 0, 1];
st{curfig}.Space=M*st{curfig}.Space;

set(handles.SIButton, 'String', state);
w_spm_orthviews('redraw');

% --- Executes on button press in PAButton.
function PAButton_Callback(hObject, eventdata, handles)
% hObject    handle to PAButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global st
curfig=handles.DPABI_fig;
curfig=w_Compatible2014bFig(curfig);
PA={'P','A'};
state=get(handles.PAButton, 'String');
state=strcmpi(state, 'P')+1;
state=PA{state};

bb=st{curfig}.bb;
offset=-bb(1,3)-bb(2,3);
M=[ 1, 0, 0, 0;
    0, -1, 0, offset;
    0, 0, 1, 0;
    0, 0, 0, 1];
st{curfig}.Space=M*st{curfig}.Space;

set(handles.PAButton, 'String', state);
w_spm_orthviews('redraw');

% --- Executes on button press in LRButton.
function LRButton_Callback(hObject, eventdata, handles)
% hObject    handle to LRButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global st
curfig=handles.DPABI_fig;
curfig=w_Compatible2014bFig(curfig);
LR={'R','L'};
state=get(handles.LRButton, 'String');
state=strcmpi(state, 'R')+1;
state=LR{state};

bb=st{curfig}.bb;
offset=-bb(1,1)-bb(2,1);
M=[-1, 0, 0, offset;
    0, 1, 0, 0;
    0, 0, 1, 0;
    0, 0, 0, 1];
st{curfig}.Space=M*st{curfig}.Space;

set(handles.LRButton, 'String', state);
w_spm_orthviews('redraw');

function index = HeaderIndex(handles, curblob)
global st;
curfig=handles.DPABI_fig;
curfig=w_Compatible2014bFig(curfig);
if nargin<2
    curblob=st{curfig}.curblob;
end

if curblob==0;
    index=0;
    return
end
OverlayFile=handles.OverlayFileName{curblob};
index=0;
for i=1:length(handles.OverlayHeaders)
    OverlayHeader=handles.OverlayHeaders{i};
    if ~isempty(OverlayHeader) && strcmp(OverlayFile, OverlayHeader.fname)
        index=i;
        break;
    end
end

if ~index
    errordlg('Header Error');
end


function centre=Pos(handles, label)
global st

centre=zeros(1,3);
switch lower(label)
    case 'xyz'
        centre(1)=str2double(get(handles.XEntry, 'String'));
        centre(2)=str2double(get(handles.YEntry, 'String'));
        centre(3)=str2double(get(handles.ZEntry, 'String'));
    case 'ijk'
        I=str2double(get(handles.IEntry, 'String'));
        J=str2double(get(handles.JEntry, 'String'));
        K=str2double(get(handles.KEntry, 'String'));
        curfig=handles.DPABI_fig;
        curfig=w_Compatible2014bFig(curfig);
        curblob=st{curfig}.curblob;
        if isfield(st{curfig}.vols{1}, 'blobs')
            tmp=st{curfig}.vols{1}.blobs{curblob}.vol.mat*[I;J;K;1];
        else
            tmp=st{curfig}.vols{1}.mat*[I;J;K;1];
        end
        centre(1)=tmp(1);
        centre(2)=tmp(2);
        centre(3)=tmp(3);
end

function OverlayHeader=SetMask(OverlayHeader)
OverlayHeader.Data(~OverlayHeader.Mask)=0;

function OverlayHeader=SetCSize(OverlayHeader)
%%
if ~OverlayHeader.CSize
    return
end
OverlayVolume=OverlayHeader.Data;
RMM=OverlayHeader.RMM;
CSize=OverlayHeader.CSize;
V=prod(OverlayHeader.Vox);
CVSize=CSize/V;

CC=bwconncomp(logical(OverlayVolume), RMM);

numPixels=cellfun(@numel, CC.PixelIdxList);

N=CC.PixelIdxList(numPixels<CVSize)';

Index=cell2mat(N);
OverlayVolume(Index)=0;

OverlayHeader.Data=OverlayVolume;

function OverlayHeader=SetPosNeg(OverlayHeader, handles)
if get(handles.OnlyPos, 'Value')
    OverlayHeader.Data(OverlayHeader.Data<0)=0;
end
if get(handles.OnlyNeg, 'Value');
    OverlayHeader.Data(OverlayHeader.Data>0)=0;
end

function OverlayVolume = SingleCluster(OverlayHeader)
OverlayVolume=OverlayHeader.Data;
RMM=OverlayHeader.RMM;

pos=w_spm_orthviews('pos');
tmp=round(inv(OverlayHeader.mat)*[pos;1]);
OI=tmp(1);
OJ=tmp(2);
OK=tmp(3);
if ~OverlayVolume(OI, OJ, OK)
    OverlayVolume=[];
    return;
end

CC=bwconncomp(logical(OverlayVolume), RMM);
L=labelmatrix(CC);
V=L(OI, OJ, OK);
OverlayVolume(L~=V)=0;

function FindPeak(OverlayHeader)
OverlayVolume=OverlayHeader.Data;
RMM=OverlayHeader.RMM;

pos=w_spm_orthviews('pos');
tmp=round(inv(OverlayHeader.mat)*[pos;1]);
OI=tmp(1);
OJ=tmp(2);
OK=tmp(3);
if ~OverlayVolume(OI, OJ, OK)
    return;
end

CC=bwconncomp(logical(OverlayVolume), RMM);
L=labelmatrix(CC);
V=L(OI, OJ, OK);
V=V==L;
OverlayVolume=OverlayVolume.*V;
Max=max(abs(OverlayVolume(:)));
[I, J, K]=ind2sub(size(OverlayVolume), find(abs(OverlayVolume)==Max));
MI=I(1);
MJ=J(1);
MK=K(1);
centre=OverlayHeader.mat*[MI;MJ;MK;1];
w_spm_orthviews('Reposition', centre(1:3));

if length(I) > 1
    msgbox(sprintf('There are %d max point in this cluster', length(I)), 'modal');
end



function AddReduce(hObject, label)
Value=str2double(get(hObject, 'String'));
if strcmpi(label, '+')
    String=num2str(Value+1);
else
    String=num2str(Value-1);
end
set(hObject, 'String', String);

% --- Executes on button press in XAddBtn.
function XAddBtn_Callback(hObject, eventdata, handles)
% hObject    handle to XAddBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AddReduce(handles.XEntry, '+');
centre=Pos(handles, 'XYZ');
w_spm_orthviews('Reposition', centre);

% --- Executes on button press in XReduceBtn.
function XReduceBtn_Callback(hObject, eventdata, handles)
% hObject    handle to XReduceBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AddReduce(handles.XEntry, '-');
centre=Pos(handles, 'XYZ');
w_spm_orthviews('Reposition', centre);

% --- Executes on button press in YAddBtn.
function YAddBtn_Callback(hObject, eventdata, handles)
% hObject    handle to YAddBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AddReduce(handles.YEntry, '+');
centre=Pos(handles, 'XYZ');
w_spm_orthviews('Reposition', centre);

% --- Executes on button press in YReduceBtn.
function YReduceBtn_Callback(hObject, eventdata, handles)
% hObject    handle to YReduceBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AddReduce(handles.YEntry, '-');
centre=Pos(handles, 'XYZ');
w_spm_orthviews('Reposition', centre);

% --- Executes on button press in ZAddBtn.
function ZAddBtn_Callback(hObject, eventdata, handles)
% hObject    handle to ZAddBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AddReduce(handles.ZEntry, '+');
centre=Pos(handles, 'XYZ');
w_spm_orthviews('Reposition', centre);

% --- Executes on button press in ZReduceBtn.
function ZReduceBtn_Callback(hObject, eventdata, handles)
% hObject    handle to ZReduceBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AddReduce(handles.ZEntry, '-');
centre=Pos(handles, 'XYZ');
w_spm_orthviews('Reposition', centre);

% --- Executes on button press in IAddBtn.
function IAddBtn_Callback(hObject, eventdata, handles)
% hObject    handle to IAddBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AddReduce(handles.IEntry, '+');
centre=Pos(handles, 'IJK');
w_spm_orthviews('Reposition', centre);

% --- Executes on button press in IReduceBtn.
function IReduceBtn_Callback(hObject, eventdata, handles)
% hObject    handle to IReduceBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AddReduce(handles.IEntry, '-');
centre=Pos(handles, 'IJK');
w_spm_orthviews('Reposition', centre);

% --- Executes on button press in JAddBtn.
function JAddBtn_Callback(hObject, eventdata, handles)
% hObject    handle to JAddBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AddReduce(handles.JEntry, '+');
centre=Pos(handles, 'IJK');
w_spm_orthviews('Reposition', centre);

% --- Executes on button press in JReduceBtn.
function JReduceBtn_Callback(hObject, eventdata, handles)
% hObject    handle to JReduceBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AddReduce(handles.JEntry, '-');
centre=Pos(handles, 'IJK');
w_spm_orthviews('Reposition', centre);

% --- Executes on button press in KAddBtn.
function KAddBtn_Callback(hObject, eventdata, handles)
% hObject    handle to KAddBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AddReduce(handles.KEntry, '+');
centre=Pos(handles, 'IJK');
w_spm_orthviews('Reposition', centre);

% --- Executes on button press in KReduceBtn.
function KReduceBtn_Callback(hObject, eventdata, handles)
% hObject    handle to KReduceBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AddReduce(handles.KEntry, '-');
centre=Pos(handles, 'IJK');
w_spm_orthviews('Reposition', centre);


function AtlasEntry_Callback(hObject, eventdata, handles)
% hObject    handle to AtlasEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AtlasEntry as text
%        str2double(get(hObject,'String')) returns contents of AtlasEntry as a double


% --- Executes during object creation, after setting all properties.
function AtlasEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AtlasEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AtlasButton.
function AtlasButton_Callback(hObject, eventdata, handles)
% hObject    handle to AtlasButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
curfig=handles.DPABI_fig;
curfig=w_Compatible2014bFig(curfig);
flag=w_AtlasSelect(curfig);
if flag
    w_spm_orthviews('Redraw');
end

% --- Executes on button press in StructrueButton.
function StructrueButton_Callback(hObject, eventdata, handles)
% hObject    handle to StructrueButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global st
curfig=handles.DPABI_fig;
curfig=w_Compatible2014bFig(curfig);
AtlasInfo=st{curfig}.AtlasInfo;
if isempty(AtlasInfo)
    msgbox('Please Select Atlas First!', 'modal');
    return
end

w_StructuralSelect(curfig);


% --- Executes during object deletion, before destroying properties.
function DPABI_fig_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to DPABI_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global st
curfig=handles.DPABI_fig;
curfig=w_Compatible2014bFig(curfig);
if ~isempty(st{curfig}.TCFlag) && ishandle(st{curfig}.TCFlag)
    delete(st{curfig}.TCFlag);
end

if ~isempty(st{curfig}.SSFlag) && ishandle(st{curfig}.SSFlag)
    delete(st{curfig}.SSFlag);
end

for i=1:numel(st{curfig}.MPFlag)
    if ~isempty(st{curfig}.MPFlag{i}) && ishandle(st{curfig}.MPFlag{i})
        delete(st{curfig}.MPFlag{i});
    end
end
st{curfig}.MPFlag=[];

st{curfig}=[];
delete(handles.DPABI_fig);


% --- Executes on button press in ModeButton.
function ModeButton_Callback(hObject, eventdata, handles)
% hObject    handle to ModeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global st
curfig=handles.DPABI_fig;
curfig=w_Compatible2014bFig(curfig);
Mode={'1','0'};
mode=get(handles.ModeButton, 'String');
state=strcmpi(mode, '1')+1;
mode=Mode{state};

set(handles.ModeButton, 'String', mode);
st{curfig}.mode=str2double(mode);
ShowUnderlay(handles);



% --- Executes on selection change in ColormapPopup.
function ColormapPopup_Callback(hObject, eventdata, handles)
% hObject    handle to ColormapPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index=HeaderIndex(handles);
OverlayHeader=handles.OverlayHeaders{index};
Value=get(handles.ColormapPopup, 'Value');
switch Value
    case 1 %JET
        ColorMap=colormap('jet(64)');
    case 2 %HSV
        ColorMap=colormap('hsv(64)');
    case 3 %HOT
        ColorMap=colormap('hot(64)');
    case 4 %Cool
        ColorMap=colormap('cool(64)');
    case 5 %Copper
        ColorMap=colormap('copper(64)');
    case 6 %Spring
        ColorMap=colormap('spring(64)');
    case 7 %Summer
        ColorMap=colormap('summer(64)');
    case 8 %Autumn
        ColorMap=colormap('autumn(64)');
    case 9 %Winter
        ColorMap=colormap('winter(64)');
    case 10%Gray
        ColorMap=colormap('gray(64)');
    case 11%Bone
        ColorMap=colormap('bone(64)');
    case 12%Pink
        ColorMap=colormap('pink(64)');
    case 13%Random 
        ColorMap=colormap('lines(64)');
    case 14%AFNI
        ColorMap=y_AFNI_ColorMap(12);
    case 15%Custom
        OverlayHeader=handles.OverlayHeaders{index};
        colormap(OverlayHeader.ColorMap);
        colormapeditor(handles.DPABI_fig);
        return;
        %ColorMap=colormap;
end

OverlayHeader.ColorMap=ColorMap;
OverlayHeader.ColorMapIndex=Value;
OverlayHeader=RedrawOverlay(OverlayHeader);
handles.OverlayHeaders{index}=OverlayHeader;
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns ColormapPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ColormapPopup


% --- Executes during object creation, after setting all properties.
function ColormapPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ColormapPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function AlphaSlider_Callback(hObject, eventdata, handles)
% hObject    handle to AlphaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AlphaListener_Callback(hObject, eventdata, hObject);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

function AlphaListener_Callback(hObject, eventdata, hFig)
handles=guidata(hFig);
curfig=handles.DPABI_fig;
curfig=w_Compatible2014bFig(curfig);

value=get(handles.AlphaSlider, 'Value');
index=HeaderIndex(handles);
if ~index
    return
end

OverlayHeader=handles.OverlayHeaders{index};
OverlayHeader.Alpha=value;
OverlayHeader=RedrawOverlay(OverlayHeader);
handles.OverlayHeaders{index}=OverlayHeader;
guidata(hFig, handles);
UpdateOverlayHandle(handles, 'On');

% --- Executes during object creation, after setting all properties.
function AlphaSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AlphaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in ColorbarTypePopup.
function ColorbarTypePopup_Callback(hObject, eventdata, handles)
% hObject    handle to ColorbarTypePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Value=get(handles.ColorbarTypePopup, 'Value');
switch Value
    case 1
        PN_Flag='+';
    case 2
        PN_Flag='-';
    case 3
        PN_Flag=[];
end
index=HeaderIndex(handles);
OverlayHeader=handles.OverlayHeaders{index};
OverlayHeader.PN_Flag=PN_Flag;
OverlayHeader=RedrawOverlay(OverlayHeader);
handles.OverlayHeaders{index}=OverlayHeader;
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns ColorbarTypePopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ColorbarTypePopup


% --- Executes during object creation, after setting all properties.
function ColorbarTypePopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ColorbarTypePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in CSTypePopup.
function CSTypePopup_Callback(hObject, eventdata, handles)
% hObject    handle to CSTypePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Value=get(handles.CSTypePopup, 'Value');
switch Value
    case 1
        RMM=6;
    case 2
        RMM=18;
    case 3
        RMM=26;
end
index=HeaderIndex(handles);
OverlayHeader=handles.OverlayHeaders{index};
OverlayHeader.RMM=RMM;
OverlayHeader=RedrawOverlay(OverlayHeader);
handles.OverlayHeaders{index}=OverlayHeader;
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns CSTypePopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CSTypePopup


% --- Executes during object creation, after setting all properties.
function CSTypePopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CSTypePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CSEntry_Callback(hObject, eventdata, handles)
% hObject    handle to CSEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CVSize=get(handles.CSEntry, 'String');
CVSize=str2double(CVSize);
index=HeaderIndex(handles);
OverlayHeader=handles.OverlayHeaders{index};

OverlayHeader.CSize=prod(OverlayHeader.Vox)*CVSize;

OverlayHeader=RedrawOverlay(OverlayHeader);
handles.OverlayHeaders{index}=OverlayHeader;
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of CSEntry as text
%        str2double(get(hObject,'String')) returns contents of CSEntry as a double


% --- Executes during object creation, after setting all properties.
function CSEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CSEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FindPeakBtn.
function FindPeakBtn_Callback(hObject, eventdata, handles)
% hObject    handle to FindPeakBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index=HeaderIndex(handles);
OverlayHeader=handles.OverlayHeaders{index};
OverlayHeader=SetMask(OverlayHeader);
OverlayHeader=SetCSize(OverlayHeader);
OverlayHeader=SetPosNeg(OverlayHeader, handles);
FindPeak(OverlayHeader);

% --- Executes on button press in CIReportBtn.
function CIReportBtn_Callback(hObject, eventdata, handles)
% hObject    handle to CIReportBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index=HeaderIndex(handles);
OverlayHeader=handles.OverlayHeaders{index};
OverlayHeader=SetMask(OverlayHeader);
OverlayHeader=SetCSize(OverlayHeader);               
OverlayHeader=SetPosNeg(OverlayHeader, handles); %Add by Sandy to set cluster size at first when someone use Cluster Report       
ClusterReport(OverlayHeader.Data, OverlayHeader, OverlayHeader.RMM);

% --- Executes on selection change in CorrectPopup.
function CorrectPopup_Callback(hObject, eventdata, handles)
% hObject    handle to CorrectPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index=HeaderIndex(handles);
if ~index
    return
end
Value=get(handles.CorrectPopup, 'Value');
OverlayHeader=handles.OverlayHeaders{index};
switch Value
    case 1
        return
    case 2 %FDR
        OverlayHeader=SeeCAT_FDRCorrection(OverlayHeader);
        if isempty(OverlayHeader)
            return
        end
        OverlayHeader=RedrawOverlay(OverlayHeader, handles.DPABI_fig);
        handles.OverlayHeaders{index}=OverlayHeader;        
    case 3 %GRF
        OverlayHeader=SeeCAT_GRFCorrection(OverlayHeader);
        if isempty(OverlayHeader)
            return
        end
        [OverlayHeader, SendHeader]=RedrawOverlay(OverlayHeader);
        %OverlayHeader.Data = SendHeader.Data; %OverlayHeader=SetCSize(OverlayHeader); %YAN Chao-Gan, 140822. Need to save the data after setting cluster size.
        handles.OverlayHeaders{index}=OverlayHeader; 
end
guidata(hObject, handles);
UpdateOverlayHandle(handles, 'On');
% Hints: contents = cellstr(get(hObject,'String')) returns CorrectPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CorrectPopup


% --- Executes during object creation, after setting all properties.
function CorrectPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CorrectPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ExportBtn.
function ExportBtn_Callback(hObject, eventdata, handles)
% hObject    handle to ExportBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Value=get(handles.ExportBtn, 'Value');
if Value==1
    return
end
switch Value 
    case 2 %Montage
        MPFlag=SeeCAT_Montage(handles.DPABI_fig);
    case 3 %Surface
        index=HeaderIndex(handles);
        OverlayHeader=handles.OverlayHeaders{index};
        OverlayHeader=SetMask(OverlayHeader);
        OverlayHeader=SetCSize(OverlayHeader);               
        OverlayHeader=SetPosNeg(OverlayHeader, handles); %Add by Sandy to set cluster size at first when someone use Cluster Report       
        CallBrainNetViewer(OverlayHeader);
end

% --- Executes on button press in SaveCCBtn.
function SaveCCBtn_Callback(hObject, eventdata, handles)
% hObject    handle to SaveCCBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index=HeaderIndex(handles);
OverlayHeader=handles.OverlayHeaders{index};
[File , Path]=uiputfile({'*.img;*.nii;*.nii.gz','Brain Image Files (*.img;*.nii;*.nii.gz)';'*.*', 'All Files (*.*)';}, ...
    'Save All Clusters' , OverlayHeader.fname);
if ~ischar(File)
    return
end
[OverlayHeader, SendHeader]=RedrawOverlay(OverlayHeader);
        
[Path, FileName, Ext]=fileparts(fullfile(Path, File));
MaskName=sprintf('%s_mask',FileName);
Data=SingleCluster(SendHeader);
if isempty(Data)
    return
end
Mask=logical(Data);
y_Write(Data, SendHeader, fullfile(Path, [FileName, Ext]));
y_Write(Mask, SendHeader, fullfile(Path, [MaskName, Ext]));


% --- Executes on button press in SaveACBtn.
function SaveACBtn_Callback(hObject, eventdata, handles)
% hObject    handle to SaveACBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index=HeaderIndex(handles);
OverlayHeader=handles.OverlayHeaders{index};
[File , Path]=uiputfile({'*.img;*.nii;*.nii.gz','Brain Image Files (*.img;*.nii;*.nii.gz)';'*.*', 'All Files (*.*)';}, ...
    'Save All Clusters' , OverlayHeader.fname);
if ~ischar(File)
    return
end
[OverlayHeader, SendHeader]=RedrawOverlay(OverlayHeader);
        
[Path, FileName, Ext]=fileparts(fullfile(Path, File));
MaskName=sprintf('%s_mask',FileName);
Data=SendHeader.Data;
Mask=logical(Data);
y_Write(Data, SendHeader, fullfile(Path, [FileName, Ext]));
y_Write(Mask, SendHeader, fullfile(Path, [MaskName, Ext]));

% --- Executes on button press in SaveAPBtn.
function SaveAPBtn_Callback(hObject, eventdata, handles)
% hObject    handle to SaveAPBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[File, Path] = uiputfile({'*.tiff';'*.jpeg';'*.png';'*.bmp'},...
    'Save Picture As');
if ~ischar(File)
    return;
end
[Path, Name, Ext]=fileparts(fullfile(Path, File));
curfig=w_Compatible2014bFig(handles.DPABI_fig);
global st
TName=sprintf('%s_Transverse', Name);
CName=sprintf('%s_Coronal', Name);
SName=sprintf('%s_Sagittal', Name);
TData=getframe(st{curfig}.vols{1}.ax{1}.ax);
CData=getframe(st{curfig}.vols{1}.ax{2}.ax);
SData=getframe(st{curfig}.vols{1}.ax{3}.ax);
 
saveas(curfig, fullfile(Path, [Name, Ext]));
%eval(['print -r300 -dtiff -noui ''',fullfile(Path, [Name, '_300dpi', Ext]),''';']); %YAN Chao-Gan, 140806.
        
imwrite(TData.cdata, fullfile(Path, [TName, Ext]));
imwrite(CData.cdata, fullfile(Path, [CName, Ext]));
imwrite(SData.cdata, fullfile(Path, [SName, Ext]));

index=HeaderIndex(handles);
OverlayHeader=handles.OverlayHeaders{index};

ColorMap=OverlayHeader.ColorMap;
L=size(ColorMap, 1);

PN_Flag=OverlayHeader.PN_Flag;
if isempty(PN_Flag)
    NMax=OverlayHeader.NMax;
    NMin=OverlayHeader.NMin;
    PMin=OverlayHeader.PMin;
    PMax=OverlayHeader.PMax;    
    if NMax==0 && NMax==NMin;
        Scale=(L:-1:L/2+1)';
    elseif PMax==0 && PMax==PMin
        Scale=(L/2:-1:1)';
    else
        Scale=(1:L)';
    end
else
    Scale=(1:L)';
end
CBName=sprintf('%s_Colorbar', Name);
imwrite(imresize(Scale, [320*L, 200]),...
    ColorMap, fullfile(Path, [CBName, Ext]));

% --- Executes on selection change in StructuralPopup.
function StructuralPopup_Callback(hObject, eventdata, handles)
% hObject    handle to StructuralPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
w_spm_orthviews('redraw', w_Compatible2014bFig(handles.DPABI_fig));
% Hints: contents = cellstr(get(hObject,'String')) returns StructuralPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from StructuralPopup


% --- Executes during object creation, after setting all properties.
function StructuralPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StructuralPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PreviewBtn.
function PreviewBtn_Callback(hObject, eventdata, handles)
% hObject    handle to PreviewBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
w_spm_orthviews('redraw', w_Compatible2014bFig(handles.DPABI_fig));
% Hint: get(hObject,'Value') returns toggle state of PreviewBtn


function AlphaEntry_Callback(hObject, eventdata, handles)
% hObject    handle to AlphaEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index=HeaderIndex(handles);
if ~index
    return
end
OverlayHeader=handles.OverlayHeaders{index};
Alpha=str2double(get(handles.AlphaEntry, 'String'));
if isnan(Alpha) || (Alpha>1 || Alpha<0) 
    errordlg('Invalid Alpha!');
    return
end
OverlayHeader.Alpha=Alpha;
handles.OverlayHeaders{index}=OverlayHeader;
guidata(hObject, handles);
UpdateOverlayHandle(handles, 'On');
% Hints: get(hObject,'String') returns contents of AlphaEntry as text
%        str2double(get(hObject,'String')) returns contents of AlphaEntry as a double


% --- Executes during object creation, after setting all properties.
function AlphaEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AlphaEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%%
% Some Plugin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function colormapeditor(fig)
%COLORMAPEDITOR starts colormap editor ui
%
%   When started the colormap editor displays the current figure's colormap
%   as a strip of rectangluar cells. Nodes, displayed as rectangles(end
%   nodes) or carrots below the strip, separate regions of uniform slope in
%   R,G and B.   The grays colormap, for example, has only two nodes since
%   R,G and B increase at a constant slope from the first to last index in
%   the map. Most of the other standard colormaps have additional nodes where
%   the slopes of the R,G or B curves change.
% 
%   As the mouse is moved over the color cells the colormap index, CData
%   value r,g,b and h,s,v values for that cell are displayed in the Current
%   Color Info box of the gui.
%
%   To add a node: 
%       Click below the cell where you wish to add the node
%
%   To delete a node:
%       Select the node by clicking on it and hit the Delete key or
%       Edit->Delete, or Ctrl-X
%
%   To move a node:
%       Click and drag or select and use left and right arrow keys.
%
%   To change a node's color:
%       Double click or click to select and Edit->Set Node Color. If multiple
%       nodes are selected the color change will apply to the last node
%       selected (current node).
%
%   To select this node, that node and all nodes between:
%       Click on this node, then Shift-Click on that node.
%
%   To select this, that and the other node:
%       Click on this, Ctrl-Click on that and the other.
%
%   To move multiple nodes at once: 
%       Select multiple nodes then use left and right arrow keys to move them
%       all at once.  Movement will stop when one of the selected nodes bumps
%       into an unselected node or end node. 
%
%   To delete multiple nodes:
%       Select the nodes and hit the Delete key, or Edit->Delete, or Ctrl-X.
%
%   To avoid flashing while editing the colormap set the figures DoubleBuffer
%   property to 'on'.
%
%   The "Interpolating Colorspace" selection determines what colorspace is
%   used to calculate the color of cells between nodes.  Initially this is
%   set to RGB, meaning that the R,G and B values of cells between nodes are
%   linearly interpolated between the R,G and B values of the nodes. Changing
%   the interpolating colorspace to HSV causes the cells between nodes to be
%   re-calculated by interpolating their H,S and V values from  the H,S and V
%   values of the nodes.  This usually produces very different results.
%   Because Hue is conceptually mapped about a color circle, the
%   interpolation between Hue values could be ambiguous.  To minimize
%   ambiguity the shortest distance around the circle is used.  For example,
%   if two  nodes have Hues of 2(slightly orange red) and 356 (slightly
%   magenta red), the cells between them would not have hues 3,4,5 ....
%   353,354,355  (orange/red-yellow-green-cyan-blue-magenta/red) but 357,
%   358, 1, 2  (orange/red-red-magenta/red).
%
%   The "Color Data Min" and "Color Data Max" editable text areas contain 
%   the values that correspond to the current axes "Clim" property.  These
%   values may be set here and are useful for selecting a range of data
%   values to which the colormap will map.  For example, your CData values
%   might range from 0 to 100, but the range of interest may be between 40 
%   and 60.  Using Color Data Min and Max (or the Axes Clim property) the
%   full variation of the colormap can be placed between the values 40 and 60
%   to improve visual/color resolution in that range.   Color Data Min
%   corresponds to Clim(0) and is the CData value to which the first Colormap
%   index is mapped.  All CData Values below this will display in the same
%   color as the first index.  Color Data Max corresponds to Clim(1) and is
%   the CData value to which the last Colormap index is mapped.  All CData 
%   values greater than this will display the same color as the last index.
%   Setting Color Data Min and Max (or Clim) will only affect the display of
%   objects (Surfaces, Patches, Images) who's CDataMapping Property is set to
%   'Scaled'.  e.g.  imagesc(im) but not image(im).
%
%   Immediate Apply is checked by default, and applies changes as they are
%   made.  If Immediate Apply is unselected, changes in the colormap editor
%   will not be applied until the apply or OK button is selected, and any
%   changes made will be lost if the Cancel button is selected.
%
%   See also COLORMAP


%   G. DeLoid 03/04/2002
%
%   Copyright 1984-2010 The MathWorks, Inc.
%   $Revision: 1.2.4.24.2.1 $  $Date: 2010/07/06 14:40:31 $

error(javachk('mwt', 'The Colormap Editor'));

import com.mathworks.page.cmapeditor.*;

% look for figure input
if nargin
    % colormapeditor([]) should do nothing
    if isempty(fig)
        return
    end
    
    if nargin == 1
        if ~any(ishandle_valid(fig)) || ~ishghandle(fig,'figure')
            error('MATLAB:colormapeditor:InvalidFigureHandle', 'Invalid figure handle');
        end
    end
else
    fig = [];
end

% get figure if not provided
if isempty(fig)
    fig = get(0,'CurrentFigure');
    if isempty(fig)
        fig=figure;
    end
end

% make sure the colormap is valid
check_colormap(get(fig,'Colormap'))



% reuse the only one if it's there
cme = get_cmapeditor();
if ~isempty(cme)
    cme.bringToFront;
    cme.setVisible;
    return;
end

% start colormapeditor and make it the only one
if feature('HGUsingMATLABClasses') % Work around G621525
    com.mathworks.page.cmapeditor.CMEditFrame.sendJavaWindowFocusEventsToMatlab(true);
end
cme = CMEditor;
set_cmapeditor(cme);
cme.init;
    
% attach the matlab callback interface so we get notified of updates
cb = handle(cme.getMatlabQue.getCallback,'callbackproperties');
set(cb,'delayedCallback',@handle_callback)

update_colormap(get(fig,'Colormap'));

ax = get(fig,'CurrentAxes');
handle_axes_change(fig,ax,false);
cme.setFigure (java(handle(fig)));
start_listeners(fig,ax);

% all set, show it now
cme.setVisible;

%----------------------------------------------------------------------%
% FUNCTIONS CALLED BY FROM JAVA EDITOR
%----------------------------------------------------------------------%
function handle_callback(callback_source, eventData) %#ok<INUSL>
% callback_source is not used
cme = get_cmapeditor();
cme_ = eventData.getEditor;
if isempty(cme) || (~isempty(cme_) && cme_.isDisposed)
    % This error was put in to catch the scenario where a callback for the cmeditor
    % is gone before a callback is handled. Since the cmeditor goes away in response
    % to a matlab event from the figure and these callbacks are coming from JAVA, there
    % is no guarantee on the order of these events. Therefore, don't error, just igonre.
    %    error('MATLAB:colormapeditor:CmeditorExpected', 'expected cmeditor to still be there')
    return;
end

if eventData.getOperation()==MLQue.CURRENT_FIGURE_UPDATE
    cme = eventData.getEditor;
end
fig = handle(cme.getFigure); % Remove java wrapper

funcode    = eventData.getOperation();
args       = eventData.getArguments();
resultSize = eventData.getResultSize;
source     = eventData.getSource;

if ~isequal(cme,cme_)
    error('MATLAB:colormapeditor:ExpectedCmeditor', 'expected cmeditor to still be there')
end

import com.mathworks.page.cmapeditor.MLQue

switch(funcode)
 case MLQue.CMAP_UPDATE
  cmap_update(fig,args);
 case MLQue.CLIM_UPDATE
  clim_update(fig,args);
 case MLQue.CMAP_STD
  source.finished(funcode, stdcmap(args, resultSize), cme);
 case MLQue.CMAP_EDITOR_OFF
  kill_listeners(fig);
 case MLQue.CHOOSE_COLOR
  source.finished(funcode, choosecolor(args), cme);
 case MLQue.CURRENT_FIGURE_UPDATE
  oldfig = handle(eventData.getEditor.getFigure);
  drawnow
  fig = get(0,'CurrentFigure');
  if ~isequal(oldfig,fig)
      oldax = [];
      if ~isempty(oldfig) && ishghandle(oldfig)
          oldax = get(oldfig,'CurrentAxes');
      end
      set_current_figure(fig,oldfig,oldax);
  end
end

%----------------------------------------------------------------------%
function cmap_update(fig,map)

%if ~ishandle_valid(fig)
if ~any(ishandle_valid(fig))
    return;
end

cmap_listen_enable(fig,'off');
set(fig,'Colormap',map);
cmap_listen_enable(fig,'on');

handles=guidata(fig);
index=SeeCAT_Viewer('HeaderIndex', handles);
OverlayHeader=handles.OverlayHeaders{index};
OverlayHeader.ColorMap=colormap;
OverlayHeader.cbarstring='0';
OverlayHeader=SeeCAT_Viewer('RedrawOverlay', OverlayHeader);

handles.OverlayHeaders{index}=OverlayHeader;
guidata(fig, handles);


%----------------------------------------------------------------------%
function clim_update(fig,lims)

%if ~ishandle_valid(fig)
if ~any(ishandle_valid(fig))
    return;
end
ax = get(fig,'CurrentAxes');
%if ~ishandle_valid(ax)
if ~any(ishandle_valid(ax))
    return;
end

cmap_listen_enable(fig,'off');
set(ax,'clim',lims);
cmap_listen_enable(fig,'on');

%----------------------------------------------------------------------%
function map=choosecolor(vals)

r = vals(1);
g = vals(2);
b = vals(3);
map=uisetcolor([r,g,b],xlate('Select Marker Color'));

%----------------------------------------------------------------------%
function map=stdcmap(maptype, mapsize)

import com.mathworks.page.cmapeditor.MLQue

switch maptype
 case MLQue.AUTUMN
  map = autumn(mapsize);
 case MLQue.BONE
  map = bone(mapsize);
 case MLQue.COLORCUBE
  map = colorcube(mapsize);
 case MLQue.COOL
  map = cool(mapsize);
 case MLQue.COPPER
  map = copper(mapsize);
 case MLQue.FLAG
  map = flag(mapsize);
 case MLQue.GRAY
  map = gray(mapsize);
 case MLQue.HOT
  map = hot(mapsize);
 case MLQue.HSV
  map = hsv(mapsize);
 case MLQue.JET
  map = jet(mapsize);
 case MLQue.LINES
  map = lines(mapsize);
 case MLQue.PINK
  map = pink(mapsize);
 case MLQue.PRISM
  map = prism(mapsize);
 case MLQue.SPRING
  map = spring(mapsize);
 case MLQue.SUMMER
  map = summer(mapsize);
 case MLQue.VGA
  map = vga;  % special case takes no size
 case MLQue.WHITE
  map = white(mapsize);
 case MLQue.WINTER
  map = winter(mapsize);
end

%----------------------------------------------------------------------%
% function cmeditor_off(fig)
% destroy any remaining listeners and remove the only one


%----------------------------------------------------------------------%
%   MATLAB listener callbacks
%----------------------------------------------------------------------%
function currentFigureChanged(hProp, eventData, oldfig, oldax) %#ok<INUSL>
% hProp is not used
if feature('HGUsingMATLABClasses')
    fig = eventData.AffectedObject.CurrentFigure;
else
    fig = eventData.NewValue;
end
get(0, 'CurrentFigure');
set_current_figure(fig,oldfig,oldax);

%----------------------------------------------------------------------%
%   Figure listener callbacks
%----------------------------------------------------------------------%
function cmapChanged(hProp, eventData, fig) %#ok<INUSL>
% hProp is not used
try
    update_colormap(get(fig,'Colormap'))
catch err
    warning(err.identifier,'%s',err.message);
end

%----------------------------------------------------------------------%
function currentAxesChanged(hProp, eventData, oldfig, oldax) %#ok<INUSL>

if feature('HGUsingMATLABClasses')
   ax = get(eventData.AffectedObject,'CurrentAxes');
else
   ax = eventData.NewValue;
end
set_current_axes(ax,oldfig,oldax);

%----------------------------------------------------------------------%
function figureDestroyed(hProp,eventData,oldfig,oldax) %#ok<INUSL>


nfigs = length(findobj(0,'type','figure','handlevisibility','on'));
% We need to check that get_cmapeditor is not empty here because when
% the test point tcolormapeditor lvlTwo_Listeners is run for
% HGUsingMATLABClasses, the call to close all closes the figure
% linked to the ColorMapEditor after the unlinked figure, so nfigs==1 %
% then this callback fires. In this case kill_listeners expects 
% that a getappdata(0,'CMEditor') is not empty, which it normally would not
% be but in the testpoint appdata(0,'CMEditor') was cleared.
if nfigs<=1 && ~isempty(get_cmapeditor)% the one being destroyed
    destroy_matlab_listeners;
    destroy_figure_listeners(oldfig);
    destroy_axes_listeners(oldax);
    kill_listeners(oldfig);
else 
    fig=get(0,'currentfigure');
    set_current_figure(fig,oldfig,oldax);
end

%----------------------------------------------------------------------%
%   Axes Listener Callbacks
%----------------------------------------------------------------------%
function climChanged(hProp, eventData, ax) %#ok<INUSL>

cme = get_cmapeditor();
if isempty(cme)
    return
end
clim = get(ax,'Clim');
cme.getModel.setColorLimits(clim,0);

%----------------------------------------------------------------------%
function axesDestroyed(hProp, eventData, oldfig, oldax) %#ok<INUSL>

cme = get_cmapeditor();
if isempty(cme)
    return;
end

fig = handle(cme.getFigure); % Remove java wrapper
if ~any(ishandle_valid(fig))
    return;
end
ax = get(fig,'currentaxes');
set_current_axes(ax,oldfig,oldax);

%----------------------------------------------------------------------%
%   Helpers
%----------------------------------------------------------------------%
function set_current_figure(fig,oldfig,oldax)

if ~any(ishandle_valid(fig)) || isequal(fig,oldfig)
    return;
end

if strncmpi (get(handle(fig),'Tag'), 'Msgbox', 6) || ...
    strcmpi (get(handle(fig),'Tag'), 'Exit') || ...
    strcmpi (get(handle(fig),'WindowStyle'), 'Modal')
    return;
end

cme = get_cmapeditor();
if isempty(cme)
    return;
end

ax = get(fig,'CurrentAxes');
% get rid of old figure listeners
destroy_figure_listeners(oldfig);
% get rid of old axes listeners
destroy_axes_listeners(oldax);
cme.setFigure (java(handle(fig)));
create_matlab_listeners(fig,ax);
update_colormap(get(fig,'Colormap'))
create_figure_listeners(fig,ax);

handle_axes_change(fig,ax,true);

%----------------------------------------------------------------------%
function set_current_axes(ax,oldfig,oldax)

if ~any(ishandle_valid(ax)) || isequal(ax,oldax)
    return;
end

fig = ancestor(ax,'figure');

% get rid of old axes listeners
destroy_axes_listeners(oldax);

% if the new axes is invalid, get out now
if ~any(ishandle_valid(ax))
    kill_listeners(oldfig);
    return;
end

create_matlab_listeners(fig,ax);
create_figure_listeners(fig,ax);

handle_axes_change(fig, ax, true);

%----------------------------------------------------------------------%
function cmap_listen_enable(fig,onoff)

% figure listeners
if ~any(ishandle_valid(fig))
    return;
end
% just cmap
if isappdata(fig,'CMEditFigListeners')
    fl = getappdata(fig,'CMEditFigListeners');
    if isobject(fl.cmapchanged)
        fl.cmapchanged.Enabled = strcmpi(onoff,'on');
    else
        set(fl.cmapchanged,'Enabled',onoff);
    end
    setappdata(fig,'CMEditFigListeners',fl);
end

% axes listeners
ax = get(fig,'CurrentAxes');
if any(ishandle_valid(ax,'CMEditAxListeners'))
    al = getappdata(ax,'CMEditAxListeners');
    if isobject(al.climchanged)
        al.climchanged.Enabled = strcmpi(onoff,'on');
    else
        set(al.climchanged,'Enabled',onoff);
    end
    setappdata(ax,'CMEditAxListeners',al);
end

%----------------------------------------------------------------------%
function start_listeners(fig,ax)

create_matlab_listeners(fig,ax);
create_figure_listeners(fig,ax);

handle_axes_change(fig,ax,true);

%----------------------------------------------------------------------%
function kill_listeners(fig)

% make sure the colormap editor is gone
cme = get_cmapeditor();
if isempty(cme)
    error('MATLAB:colormapeditor:ColormapeditorAppdataExpected',...
        'expected colormapeditor appdata to still be there')
end
cme.close;

% we need to kill these now, otherwise we'll leak the listeners and
% they will continue to fire after this colormap editor is gone
destroy_matlab_listeners

if any(ishandle_valid(fig))
    destroy_figure_listeners(fig);

    % axes
    ax = get(fig,'CurrentAxes');

    % return if no current axes or it is being destroyed
    if any(ishandle_valid(ax))
        destroy_axes_listeners(ax);
    end
end

% now flush out the cmap editor handle
rm_cmapeditor();

%----------------------------------------------------------------------%
function create_matlab_listeners(fig,ax)

if feature('HGUsingMATLABClasses')
   rt = handle(0);
   ml.cfigchanged = event.proplistener(rt,rt.findprop('CurrentFigure'), ...
        'PostSet',@(es,ed) currentFigureChanged(es,ed,fig,ax));
else
    cls = classhandle(handle(0));
    ml.cfigchanged = handle.listener(0, cls.findprop('CurrentFigure'), ...
        'PropertyPostSet', {@currentFigureChanged, fig, ax});
end
setappdata(0,'CMEditMATLABListeners',ml);

%----------------------------------------------------------------------%
function destroy_matlab_listeners


if isappdata(0,'CMEditMATLABListeners');
    % we actually need to delete these handles or they
    % will continue to fire
    ld = getappdata(0,'CMEditMATLABListeners');
    fn = fields(ld);
    for i = 1:length(fn)
        l = ld.(fn{i});
        if ishghandle(l)
            delete(l);
        end
    end
    rmappdata(0,'CMEditMATLABListeners');
end

%----------------------------------------------------------------------%
function create_figure_listeners(fig,ax)

if any(ishandle_valid(fig))
    
    fig = handle(fig);
    if feature('HGUsingMATLABClasses')
        fl.deleting = event.listener(fig, ...
                  'ObjectBeingDestroyed', @(es,ed) figureDestroyed(es,ed,fig, ax));
        fl.cmapchanged = event.proplistener(fig,fig.findprop('Colormap'), ...
                  'PostSet',@(es,ed) cmapChanged(es,ed,fig));
        fl.caxchanged = event.proplistener(fig, fig.findprop('CurrentAxes'), ...
                  'PostSet',@(es,ed) currentAxesChanged(es,ed,fig,ax));
    else     
        cls = classhandle(handle(fig));
        fl.deleting = handle.listener(fig, ...
                  'ObjectBeingDestroyed', {@figureDestroyed,fig, ax});
        fl.cmapchanged = handle.listener(fig, cls.findprop('Colormap'), ...
                  'PropertyPostSet', {@cmapChanged, fig});
        fl.caxchanged = handle.listener(fig, cls.findprop('CurrentAxes'), ...
                  'PropertyPostSet', {@currentAxesChanged, fig, ax});
    end
    setappdata(fig,'CMEditFigListeners',fl);
end

%----------------------------------------------------------------------%
function enable_figure_listeners(fig,onoff)

if any(ishandle_valid(fig, 'CMEditFigListeners'))
    fl = getappdata(fig,'CMEditFigListeners');
    if isobject(fl.cmapchanged)
        fl.cmapchanged.Enabled = strcmpi(onoff,'on');
    else
        set(fl.cmapchanged,'Enabled',onoff);
    end
    if isobject(fl.caxchanged)
        fl.caxchanged.Enabled = strcmpi(onoff,'on');
    else
        set(fl.caxchanged,'Enabled',onoff);
    end
    if isobject(fl.deleting)
        fl.deleting.Enabled = strcmpi(onoff,'on');
    else
        set(fl.deleting,'Enabled',onoff);
    end
    setappdata(fig,'CMEditFigListeners',fl);
end

%----------------------------------------------------------------------%
function destroy_figure_listeners(fig)

enable_figure_listeners(fig,'off');
if any(ishandle_valid(fig, 'CMEditFigListeners'))
    rmappdata(fig,'CMEditFigListeners');
end

%----------------------------------------------------------------------%
function create_axes_listeners(fig,ax)

if any(ishandle_valid(ax))
    if feature('HGUsingMATLABClasses')
        al.deleting = event.listener(ax, ...
                  'ObjectBeingDestroyed',@(es,ed) axesDestroyed(es,ed,fig,ax));
        al.climchanged = event.proplistener(ax,ax.findprop('Clim'), ...
                  'PostSet', @(es,ed) climChanged(es,ed,ax));
    else
        cls = classhandle(handle(ax));
        ax = handle(ax);
        al.deleting = handle.listener(ax, ...
                  'ObjectBeingDestroyed',{@axesDestroyed,fig,ax});
        al.climchanged = handle.listener(ax, cls.findprop('Clim'), ...
                  'PropertyPostSet', {@climChanged, ax});
    end
    setappdata(ax,'CMEditAxListeners',al);
end

%----------------------------------------------------------------------%
function enable_axes_listeners(ax,onoff)

if any(ishandle_valid(ax, 'CMEditAxListeners'))
    al = getappdata(ax,'CMEditAxListeners');
    if isobject(al.climchanged)
        al.climchanged.Enabled = strcmpi(onoff,'on');
    else
        set(al.climchanged,'Enabled',onoff);
    end
    if isobject(al.deleting)
        al.deleting.Enabled = strcmpi(onoff,'on');
    else
        set(al.deleting,'Enabled',onoff);
    end
    setappdata(ax,'CMEditAxListeners',al);
end

%----------------------------------------------------------------------%
function destroy_axes_listeners(ax)

enable_axes_listeners(ax,'off');
if any(ishandle_valid(ax, 'CMEditAxListeners'))
    rmappdata(ax,'CMEditAxListeners');
end

%----------------------------------------------------------------------%
function update_colormap(cmap)
check_colormap(cmap);
cme = get_cmapeditor();
if ~isempty(cme) && ~isempty(cme.getModel) 
    % cme.getModel.setColorMapModel(cmap);
    cme.getModel.setBestColorMapModel(cmap);
end

%----------------------------------------------------------------------%
function yesno = ishandle_valid(h,appdata_field)
error(nargchk(1,2,nargin,'struct'));

if nargin == 1
    appdata_field = [];
end
yesno = any(ishghandle(h)) && ~strcmpi('on',get(h,'BeingDeleted'));
if yesno && ~isempty(appdata_field)
    yesno = yesno && isappdata(h,appdata_field);
end        

%----------------------------------------------------------------------%
function handle_axes_change(fig,ax,create_listeners)
cme = get_cmapeditor();
if isempty(cme)
    return;
end
if isempty(cme.getFrame) || isempty(cme.getModel)
    return;
end

if ~any(ishandle_valid(ax))
    cme.getFrame.setColorLimitsEnabled(0);
else
    clim = get(ax,'Clim');
    cme.getFrame.setColorLimitsEnabled(1);
    cme.getModel.setColorLimits(clim,0);
    if (create_listeners)
        create_axes_listeners(fig,ax);
    end
end

%----------------------------------------------------------------------%
function check_colormap(cmap)
if isempty(cmap)
    error('MATLAB:colormapeditor:ColormapEmpty', 'Empty colormap');
end

%----------------------------------------------------------------------%
function cme = get_cmapeditor
cme = getappdata(0,'CMEditor');

%----------------------------------------------------------------------%
function set_cmapeditor(cme)
setappdata(0,'CMEditor',cme);

%----------------------------------------------------------------------%
function rm_cmapeditor
rmappdata(0,'CMEditor');

function ClusterReport(data,head,ClusterConnectivityCriterion)
% Generate the report of the thresholded clusters.
% Based on CUI Xu's xjview. (http://www.alivelearn.net/xjview/)
% Revised by YAN Chao-Gan and ZHU Wei-Xuan 20091108: suitable for different Cluster Connectivity Criterion: surface connected, edge connected, corner connected.

if ~(exist('TDdatabase.mat'))
    uiwait(msgbox('This function is based on CUI Xu''s xjview, please install xjview8 or later version at first (http://www.alivelearn.net/xjview/).','REST Slice Viewer'));
    return
end

disp('This report is based on CUI Xu''s xjview. (http://www.alivelearn.net/xjview/)'); 
disp('Revised by YAN Chao-Gan and ZHU Wei-Xuan 20091108: suitable for different Cluster Connectivity Criterion: surface connected, edge connected, corner connected.');
nozeropos=find(data~=0);
[i j k]=ind2sub(size(data),nozeropos);
cor=[i j k];
mni=cor2mni(cor,head.mat);

if isempty(mni)
    %errordlg('No cluster is picked up.','oops');
    disp( 'No cluster is found. So no report will be generated.'); 
    return;
end

intensity=data(nozeropos);


L=cor';
dim = [max(L(1,:)) max(L(2,:)) max(L(3,:))];
vol = zeros(dim(1),dim(2),dim(3));
indx = sub2ind(dim,L(1,:)',L(2,:)',L(3,:)');
vol(indx) = 1;
[cci,num] = bwlabeln(vol,ClusterConnectivityCriterion);
A = cci(indx');

clusterID = unique(A);
numClusters = length(clusterID);
disp(['Number of clusters found: ' num2str(numClusters)]);

h=waitbar(0, 'Please wait...');
for mm = clusterID
    pos = find(A == clusterID(mm));
    numVoxels = length(pos);
    tmpmni = mni(pos,:);
    tmpintensity = intensity(pos);
    
    peakpos = find(abs(tmpintensity) == max(abs(tmpintensity)));
    peakcoord = tmpmni(peakpos,:);
    peakintensity = tmpintensity(peakpos);
    
        % list structure of voxels in this cluster
    x = load('TDdatabase.mat');
    [a, b] = cuixuFindStructure(tmpmni, x.DB);
    names = unique(b(:));
    index = NaN*zeros(length(b(:)),1);
    for ii=1:length(names)
        pos = find(strcmp(b(:),names{ii}));
        index(pos) = ii;
    end

    report = {};
    
    for ii=1:max(index)
        report{ii,1} = names{ii};
        report{ii,2} = length(find(index==ii));
    end
    for ii=1:size(report,1)
        for jj=ii+1:size(report,1)
            if report{ii,2} < report{jj,2}
                tmp = report(ii,:);
                report(ii,:) = report(jj,:);
                report(jj,:) = tmp;
            end
        end
    end
    report = [{'structure','# voxels'}; {'--TOTAL # VOXELS--', length(a)}; report];

    report2 = {sprintf('%s\t%s',report{1,2}, report{1,1}),''};
    for ii=2:size(report,1)
        if strcmp('undefined', report{ii,1}); continue; end
        report2 = [report2, {sprintf('%5d\t%s',report{ii,2}, report{ii,1})}];
    end

    disp(['----------------------'])
    disp(['Cluster ' num2str(mm)])
    disp(['Number of voxels: ' num2str(numVoxels)])
    
    if size(peakcoord,1)<=1; %YAN Chao-Gan, 100814. If multi-voxels have the same peak value, then skip display the peak information.
        disp(['Peak MNI coordinate: ' num2str(peakcoord)])
        [a,b] = cuixuFindStructure(peakcoord, x.DB);
        disp(['Peak MNI coordinate region: ' a{1}]);
        disp(['Peak intensity: ' num2str(peakintensity)])
    end
    
    for kk=1:length(report2)
        disp(report2{kk});
    end
    waitbar(mm/clusterID(end), h);
end
close(h);
return

function mni = cor2mni(cor, T)
% function mni = cor2mni(cor, T)
% convert matrix coordinate to mni coordinate
%
% cor: an Nx3 matrix
% T: (optional) rotation matrix
% mni is the returned coordinate in mni space
%
% caution: if T is not given, the default T is
% T = ...
%     [-4     0     0    84;...
%      0     4     0  -116;...
%      0     0     4   -56;...
%      0     0     0     1];
%
% xu cui
% 2004-8-18
% last revised: 2005-04-30

if nargin == 1
    T = ...
        [-4     0     0    84;...
         0     4     0  -116;...
         0     0     4   -56;...
         0     0     0     1];
end

cor = round(cor);
mni = T*[cor(:,1) cor(:,2) cor(:,3) ones(size(cor,1),1)]';
mni = mni';
mni(:,4) = [];
return;

function coordinate = mni2cor(mni, T)
% function coordinate = mni2cor(mni, T)
% convert mni coordinate to matrix coordinate
%
% mni: a Nx3 matrix of mni coordinate
% T: (optional) transform matrix
% coordinate is the returned coordinate in matrix
%
% caution: if T is not specified, we use:
% T = ...
%     [-4     0     0    84;...
%      0     4     0  -116;...
%      0     0     4   -56;...
%      0     0     0     1];
%
% xu cui
% 2004-8-18
%

if isempty(mni)
    coordinate = [];
    return;
end

if nargin == 1
	T = ...
        [-4     0     0    84;...
         0     4     0  -116;...
         0     0     4   -56;...
         0     0     0     1];
end

coordinate = [mni(:,1) mni(:,2) mni(:,3) ones(size(mni,1),1)]*(inv(T))';
coordinate(:,4) = [];
coordinate = round(coordinate);
return;

function [onelinestructure, cellarraystructure] = cuixuFindStructure(mni, DB)
% function [onelinestructure, cellarraystructure] = cuixuFindStructure(mni, DB)
%
% this function converts MNI coordinate to a description of brain structure
% in aal
%
%   mni: the coordinates (MNI) of some points, in mm.  It is Nx3 matrix
%   where each row is the coordinate for one point
%   LDB: the database.  This variable is available if you load
%   TDdatabase.mat
%
%   onelinestructure: description of the position, one line for each point
%   cellarraystructure: description of the position, a cell array for each point
%
%   Example:
%   cuixuFindStructure([72 -34 -2; 50 22 0], DB)
%
% Xu Cui
% 2007-11-20
%

N = size(mni, 1);

% round the coordinates
mni = round(mni/2) * 2;

T = [...
     2     0     0   -92
     0     2     0  -128
     0     0     2   -74
     0     0     0     1];

index = mni2cor(mni, T);

cellarraystructure = cell(N, length(DB));
onelinestructure = cell(N, 1);

for ii=1:N
    for jj=1:length(DB)
        graylevel = DB{jj}.mnilist(index(ii, 1), index(ii, 2),index(ii, 3));
        if graylevel == 0
            thelabel = 'undefined';
        else
            if jj==length(DB); tmp = ' (aal)'; else tmp = ''; end
            thelabel = [DB{jj}.anatomy{graylevel} tmp];
        end
        cellarraystructure{ii, jj} = thelabel;
        onelinestructure{ii} = [ onelinestructure{ii} ' // ' thelabel ];
    end
end

function NewColorMap = AdjustColorMap(OriginalColorMap,NullColor,NMax,NMin,PMin,PMax, PN_Flag)
% Adjust the colormap to leave blank to values under threshold, the orginal color map with be set into [NMax NMin] and [PMin PMax].
% Input: OriginalColorMap - the original color map
%        NullColor - The values between NMin and PMin will be set to this color (leave blank)
%        NMax, NMin, PMin, PMax - set the axis of colorbar (the orginal color map with be set into [NMax NMin] and [PMin PMax])
% Output: NewColorMap - the generated color map, a 100000 by 3 matrix.
%___________________________________________________________________________
% Written by YAN Chao-Gan 111023.
% The Nathan Kline Institute for Psychiatric Research, 140 Old Orangeburg Road, Orangeburg, NY 10962, USA
% Child Mind Institute, 445 Park Avenue, New York, NY 10022, USA
% The Phyllis Green and Randolph Cowen Institute for Pediatric Neuroscience, New York University Child Study Center, New York, NY 10016, USA
% ycg.yan@gmail.com

NewColorMap = repmat(NullColor,[100000 1]);
ColorLen=size(OriginalColorMap,1);

%Add Full Colormap for Positive or Negative Value  
if exist('PN_Flag', 'var')==1
    if ~isempty(PN_Flag)
        if strcmpi(PN_Flag, '+')
            if PMax==PMin
                PMin=PMin-realmin;
            end
            PositiveColorSegment = ceil(100000*(PMax-PMin)/(PMax-NMax)/(ColorLen));
            for iColor=ColorLen:-1:1
                Segment=PositiveColorSegment;
                Begin=100000-(ColorLen-iColor+1)*PositiveColorSegment+1;
                if Begin < 1
                    Segment=Segment-(1-Begin);
                    Begin=1;        
                end
                End=100000-(ColorLen-iColor)*PositiveColorSegment;
                if End > 100000
                    Segment=Segment-(End-100000);
                    End=100000;
                end
                NewColorMap(Begin:End,:) = repmat(OriginalColorMap(iColor,:),[Segment 1]);
            end            
        elseif strcmpi(PN_Flag, '-')
            NegativeColorSegment = ceil(100000*(NMin-NMax)/(PMax-NMax)/(ColorLen));
            for iColor=1:ColorLen
                Segment=NegativeColorSegment;
                Begin=(iColor-1)*NegativeColorSegment+1;
                if Begin < 1
                    Segment=Segment-(1-Begin);
                    Begin=1;
                end
                End=(iColor)*NegativeColorSegment;
                if End > 100000
                    Segment=Segment-(End-100000);
                    End=100000;
                end
                NewColorMap(Begin:End,:) = repmat(OriginalColorMap(iColor,:),[Segment 1]);
            end
        end
        
        return
    end
end

NegativeColorSegment = ceil(100000*(NMin-NMax)/(PMax-NMax)/(ColorLen/2));
for iColor=1:fix(ColorLen/2)
    Segment=NegativeColorSegment;
    Begin=(iColor-1)*NegativeColorSegment+1;
    if Begin < 1
        Segment=Segment-(1-Begin);
        Begin=1;
    end
    End=(iColor)*NegativeColorSegment;
    if End > 100000
        Segment=Segment-(End-100000);
        End=100000;
    end
    NewColorMap(Begin:End,:) = repmat(OriginalColorMap(iColor,:),[Segment 1]);
end

if PMax==PMin
    PMin=PMin-realmin;
end
PositiveColorSegment = ceil(100000*(PMax-PMin)/(PMax-NMax)/(ColorLen/2));
for iColor=ColorLen:-1:ceil(ColorLen/2+1)
    Segment=PositiveColorSegment;
    Begin=100000-(ColorLen-iColor+1)*PositiveColorSegment+1;
    if Begin < 1
        Segment=Segment-(1-Begin);
        Begin=1;        
    end
    End=100000-(ColorLen-iColor)*PositiveColorSegment;
    if End > 100000
        Segment=Segment-(End-100000);
        End=100000;
    end
    NewColorMap(Begin:End,:) = repmat(OriginalColorMap(iColor,:),[Segment 1]);
end

function CallBrainNetViewer(OverlayHeader, SurfFileName, ViewType)
%% Revised from y_CallBrainNetViewer in DPABI, by Sandy Wang
[H_BrainNet] = BrainNet;
global FLAG
global EC
global surf

% Reading Surf Data. Referenced from Mingrui Xia's BrainNet Viewer
if ~exist('SurfFileName','var')
    [BrainNetViewerPath, fileN, extn] = fileparts(which('BrainNet.m'));
    SurfFileName=[BrainNetViewerPath,filesep,'Data',filesep,'SurfTemplate',filesep,'BrainMesh_ICBM152_smoothed.nv'];
end
fid=fopen(SurfFileName);
surf.vertex_number=fscanf(fid,'%f',1);
surf.coord=fscanf(fid,'%f',[3,surf.vertex_number]);
surf.ntri=fscanf(fid,'%f',1);
surf.tri=fscanf(fid,'%d',[3,surf.ntri])';
fclose(fid);

% Set up View type
if ~exist('viewtype','var')
    ViewType='MediumView';
end
if strcmpi(ViewType,'FullView')
    EC.lot.view=2;
elseif strcmpi(ViewType,'MediumView')
    EC.lot.view=3;
elseif strcmpi(ViewType,'SagittalView')
    EC.lot.view=1;
    EC.lot.view_direction=1;
elseif strcmpi(ViewType,'AxialView')
    EC.lot.view=1;
    EC.lot.view_direction=2;
elseif strcmpi(ViewType,'CoronalView')
    EC.lot.view=1;
    EC.lot.view_direction=3;
end

ColorMap=OverlayHeader.ColorMap;
NMax=OverlayHeader.NMax;
NMin=OverlayHeader.NMin;
PMin=OverlayHeader.PMin;
PMax=OverlayHeader.PMax;
if NMax==0
    NMin=-0.1;
    NMax=-10;
end
if PMax==0;
    PMin=0.1;
    PMax=10;
end

PN_Flag=OverlayHeader.PN_Flag;

% Surface Data
surf.hdr=OverlayHeader;
surf.mask=OverlayHeader.Data;
surf.vol=OverlayHeader.Data; 

% ColorMap
EC.vol.px=PMax;
EC.vol.nx=NMax;
EC.vol.pn=PMin;
EC.vol.nn=NMin;
EC.vol.CM=ColorMap;
EC.vol.CMt=ColorMap;
EC.vol.color_map=24;
EC.vol.adjustCM=1;
%EC.vol.CM=y_AdjustColorMap(ColorMap,EC.vol.null,NMax,NMin,PMin,PMax, PN_Flag);
if isempty(PN_Flag)
   EC.vol.display=1;
elseif strcmpi(PN_Flag, '+')
   EC.vol.display=2;
else
   EC.vol.display=3;
end
EC.vol.mapalgorithm=5;

% Set up other parameters
EC.msh.alpha=1;
FLAG.MAP=2;
FLAG.LF=1;
FLAG.Loadfile=9;

% Tell Brain Net Viewer is called by REST and do not reset colormap
%FLAG.EC=1;
%FLAG.IsCalledByREST=1;

% Call Brain Net Viewer ReDraw callback to refresh

set(H_BrainNet,'handlevisib','on');
BrainNet('NV_m_nm_Callback',H_BrainNet, [], guidata(H_BrainNet));

function varargout = w_spm_orthviews(action,varargin)
% Display orthogonal views of a set of images
% FORMAT H = w_spm_orthviews('Image',filename[,position])
% filename - name of image to display
% area     - position of image {relative}
%            [left, bottom, width, height]
% H        - handle for orthogonal sections
%
% FORMAT w_spm_orthviews('Reposition',centre)
% centre   - X, Y & Z coordinates of centre voxel
%
% FORMAT w_spm_orthviews('Space'[,handle[,M,dim]])
% handle   - the view to define the space by, optionally with extra
%            transformation matrix and dimensions (e.g. one of the blobs
%            of a view)
% with no arguments - puts things into mm space
%
% FORMAT H = w_spm_orthviews('Caption', handle, string, [Property, Value])
% handle   - the view to which a caption should be added
% string   - the caption text to add
% optional:  Property-Value pairs, e.g. 'FontWeight', 'Bold'
%
% H        - the handle to the object whose String property has the caption
%
% FORMAT w_spm_orthviews('BB',bb)
% bb       - bounding box
%            [loX loY loZ
%             hiX hiY hiZ]
%
% FORMAT w_spm_orthviews('MaxBB')
% Set the bounding box big enough to display the whole of all images
%
% FORMAT w_spm_orthviews('Resolution'[,res])
% res      - resolution (mm)
% Set the sampling resolution for all images. The effective resolution
% will be the minimum of res and the voxel sizes of all images. If no
% resolution is specified, the minimum of 1mm and the voxel sizes of the
% images is used.
%
% FORMAT w_spm_orthviews('Zoom'[,fov[,res]])
% fov      - half width of field of view (mm)
% res      - resolution (mm)
% Set the displayed part and sampling resolution for all images. The
% image display will be centered at the current crosshair position. The
% image region [xhairs-fov xhairs+fov] will be shown.
% If no argument is given or fov == Inf, the image display will be reset to
% "Full Volume". If fov == 0, the image will be zoomed to the bounding box
% from spm_get_bbox for the non-zero voxels of the image. If fov is NaN,
% then a threshold can be entered, and spm_get_bbox will be used to derive
% the bounding box of the voxels above this threshold.
% Optionally, the display resolution can be set as well.
%
% FORMAT w_spm_orthviews('Redraw')
% Redraw the images
%
% FORMAT w_spm_orthviews('Reload_mats')
% Reload the voxel-world mapping matrices from the headers stored on disk,
% e.g. following reorientation of some images.
%
% FORMAT w_spm_orthviews('Delete', handle)
% handle   - image number to delete
%
% FORMAT w_spm_orthviews('Reset')
% Clear the orthogonal views
%
% FORMAT w_spm_orthviews('Pos')
% Return the co-ordinate of the crosshairs in millimetres in the
% standard space.
%
% FORMAT w_spm_orthviews('Pos', i)
% Return the voxel co-ordinate of the crosshairs in the image in the
% ith orthogonal section.
%
% FORMAT w_spm_orthviews('Xhairs','off') OR w_spm_orthviews('Xhairs')
% Disable the cross-hairs on the display
%
% FORMAT w_spm_orthviews('Xhairs','on')
% Enable the cross-hairs
%
% FORMAT w_spm_orthviews('Interp',hld)
% Set the hold value to hld (see spm_slice_vol)
%
% FORMAT w_spm_orthviews('AddBlobs',handle,XYZ,Z,mat,name)
% Add blobs from a pointlist to the image specified by the handle(s)
% handle   - image number to add blobs to
% XYZ      - blob voxel locations
% Z        - blob voxel intensities
% mat      - matrix from voxels to millimeters of blob.
% name     - a name for this blob
% This method only adds one set of blobs, and displays them using a split
% colour table.
%
% FORMAT w_spm_orthviews('SetBlobsMax', vn, bn, mx)
% Set maximum value for blobs overlay number bn of view number vn to mx.
%
% FORMAT w_spm_orthviews('AddColouredBlobs',handle,XYZ,Z,mat,colour,name)
% Add blobs from a pointlist to the image specified by the handle(s)
% handle   - image number to add blobs to
% XYZ      - blob voxel locations
% Z        - blob voxel intensities
% mat      - matrix from voxels to millimeters of blob.
% colour   - the 3 vector containing the colour that the blobs should be
% name     - a name for this blob
% Several sets of blobs can be added in this way, and it uses full colour.
% Although it may not be particularly attractive on the screen, the colour
% blobs print well.
%
% FORMAT w_spm_orthviews('AddColourBar',handle,blobno)
% Add colourbar for a specified blob set
% handle    - image number
% blobno    - blob number
%
% FORMAT w_spm_orthviews('RemoveBlobs',handle)
% Remove all blobs from the image specified by the handle(s)
%
% FORMAT w_spm_orthviews('Addtruecolourimage',handle,filename,colourmap,prop,mx,mn)
% Add blobs from an image in true colour
% handle    - image number to add blobs to [Default: 1]
% filename  - image containing blob data [Default: GUI input]
% colourmap - colormap to display blobs in [Default: GUI input]
% prop      - intensity proportion of activation cf grayscale [default: 0.4]
% mx        - maximum intensity to scale to [maximum value in activation image]
% mn        - minimum intensity to scale to [minimum value in activation image]
%
% FORMAT w_spm_orthviews('Register',hReg)
% hReg      - Handle of HandleGraphics object to build registry in
% See spm_XYZreg for more information.
%
% FORMAT w_spm_orthviews('AddContext',handle)
% FORMAT w_spm_orthviews('RemoveContext',handle)
% handle    - image number to add/remove context menu to
%
% FORMAT w_spm_orthviews('ZoomMenu',zoom,res)
% FORMAT [zoom, res] = w_spm_orthviews('ZoomMenu')
% zoom      - A list of predefined zoom values
% res       - A list of predefined resolutions
% This list is used by spm_image and w_spm_orthviews('addcontext',...) to
% create the 'Zoom' menu. The values can be retrieved by calling
% w_spm_orthviews('ZoomMenu') with 2 output arguments. Values of 0, NaN and
% Inf are treated specially, see the help for w_spm_orthviews('Zoom' ...).
%__________________________________________________________________________
%
% PLUGINS
% The display capabilities of w_spm_orthviews can be extended with plugins.
% These are located in the w_spm_orthviews subdirectory of the SPM
% distribution.
% The functionality of plugins can be accessed via calls to
% w_spm_orthviews('plugin_name', plugin_arguments). For detailed descriptions
% of each plugin see help w_spm_orthviews/spm_ov_'plugin_name'.
%__________________________________________________________________________
% Copyright (C) 1996-2012 Wellcome Trust Centre for Neuroimaging

% John Ashburner et al
% $Id: w_spm_orthviews.m 5450 2013-04-26 11:25:36Z guillaume $


% The basic fields of st are:
%         n        - the number of images currently being displayed
%         vols     - a cell array containing the data on each of the
%                    displayed images.
%         Space    - a mapping between the displayed images and the
%                    mm space of each image.
%         bb       - the bounding box of the displayed images.
%         centre   - the current centre of the orthogonal views
%         callback - a callback to be evaluated on a button-click.
%         xhairs   - crosshairs off/on
%         hld      - the interpolation method
%         fig      - the figure that everything is displayed in
%         mode     - the position/orientation of the sagittal view.
%                    - currently always 1
%
%         st{curfig}.registry.hReg \_ See spm_XYZreg for documentation
%         st{curfig}.registry.hMe  /
%
% For each of the displayed images, there is a non-empty entry in the
% vols cell array.  Handles returned by "w_spm_orthviews('Image',.....)"
% indicate the position in the cell array of the newly created ortho-view.
% Operations on each ortho-view require the handle to be passed.
%
% When a new image is displayed, the cell entry contains the information
% returned by spm_vol (type help spm_vol for more info).  In addition,
% there are a few other fields, some of which are documented here:
%
%         premul  - a matrix to premultiply the .mat field by.  Useful
%                   for re-orienting images.
%         window  - either 'auto' or an intensity range to display the
%                   image with.
%         mapping - Mapping of image intensities to grey values. Currently
%                   one of 'linear', 'histeq', loghisteq',
%                   'quadhisteq'. Default is 'linear'.
%                   Histogram equalisation depends on the image toolbox
%                   and is only available if there is a license available
%                   for it.
%         ax      - a cell array containing an element for the three
%                   views.  The fields of each element are handles for
%                   the axis, image and crosshairs.
%
%         blobs   - optional.  Is there for using to superimpose blobs.
%                   vol     - 3D array of image data
%                   mat     - a mapping from vox-to-mm (see spm_vol, or
%                             help on image formats).
%                   max     - maximum intensity for scaling to.  If it
%                             does not exist, then images are auto-scaled.
%
%                   There are two colouring modes: full colour, and split
%                   colour.  When using full colour, there should be a
%                   'colour' field for each cell element.  When using
%                   split colourscale, there is a handle for the colorbar
%                   axis.
%
%                   colour  - if it exists it contains the
%                             red,green,blue that the blobs should be
%                             displayed in.
%                   cbar    - handle for colorbar (for split colourscale).
%
% PLUGINS
% The plugin concept has been developed to extend the display capabilities
% of w_spm_orthviews without the need to rewrite parts of it. Interaction
% between w_spm_orthviews and plugins takes place
% a) at startup: The subfunction 'reset_st' looks for folders
%                'w_spm_orthviews' in spm('Dir') and each toolbox
%                folder. Files with a name spm_ov_PLUGINNAME.m in any of
%                these folders will be treated as plugins.
%                For each such file, PLUGINNAME will be added to the list
%                st{curfig}.plugins{:}.
%                The subfunction 'add_context' calls each plugin with
%                feval(['spm_ov_', st{curfig}.plugins{k}], ...
%                  'context_menu', i, parent_menu)
%                Each plugin may add its own submenu to the context
%                menu.
% b) at redraw:  After images and blobs of st{curfig}.vols{i} are drawn, the
%                struct st{curfig}.vols{i} is checked for field names that occur in
%                the plugin list st{curfig}.plugins{:}. For each matching entry, the
%                corresponding plugin is called with the command 'redraw':
%                feval(['spm_ov_', st{curfig}.plugins{k}], ...
%                  'redraw', i, TM0, TD, CM0, CD, SM0, SD);
%                The values of TM0, TD, CM0, CD, SM0, SD are defined in the
%                same way as in the redraw subfunction of w_spm_orthviews.
%                It is up to the plugin to do all necessary redraw
%                operations for its display contents. Each displayed item
%                must have set its property 'HitTest' to 'off' to let events
%                go through to the underlying axis, which is responsible for
%                callback handling. The order in which plugins are called is
%                undefined.
%
%
%
%___________________________________________________________________________
% Revised by YAN Chao-Gan, 130609. 
% The Nathan Kline Institute for Psychiatric Research, 140 Old Orangeburg Road, Orangeburg, NY 10962, USA
% Child Mind Institute, 445 Park Avenue, New York, NY 10022, USA
% The Phyllis Green and Randolph Cowen Institute for Pediatric Neuroscience, New York University Child Study Center, New York, NY 10016, USA
% ycg.yan@gmail.com



global st

if nargin>2
    if ishandle(varargin{1})
        curfig=varargin{1};
    elseif all(ishandle(varargin{2}))
        curfig=varargin{2};
    end
else
    curfig=GetCurFig;
end

persistent zoomlist reslist

if isempty(st), reset_st; end

if ~nargin, action = ''; end

if ~any(strcmpi(action,{'reposition','pos'}))
    spm('Pointer','Watch');
end
    
switch lower(action)
    case 'image'
        H = specify_image(varargin{1});
        if ~isempty(H)
            if numel(varargin)>=2
                st{curfig}.vols{H}.area = varargin{2};
            else
                st{curfig}.vols{H}.area = [0 0 1 1];
            end
            if isempty(st{curfig}.bb)
                st{curfig}.bb = maxbb;
            else
                bb=maxbb;
                if ~isempty(find(st{curfig}.bb-bb, 1))
                    st{curfig}.bb=bb;
                    mmcentre     = mean(st{curfig}.Space*[maxbb';1 1],2)';
                    st{curfig}.centre    = mmcentre(1:3);
                end
            end
            resolution;
            bbox;
            cm_pos;
        end
        varargout{1} = curfig;
        if ~isfield(st{curfig}, 'centre')
            mmcentre     = mean(st{curfig}.Space*[maxbb';1 1],2)';
            st{curfig}.centre    = mmcentre(1:3);
        end
        redraw_all

    case 'caption'
        vh = valid_handles(varargin{1});
        nh = numel(vh);
        
        xlh = nan(nh, 1);
        for i = 1:nh
            xlh(i) = get(st{curfig}.vols{vh(i)}.ax{3}.ax, 'XLabel');
            if iscell(varargin{2})
                if i <= length(varargin{2})
                    set(xlh(i), 'String', varargin{2}{i});
                end
            else
                set(xlh(i), 'String', varargin{2});
            end
            for np = 4:2:nargin
                property = varargin{np-1};
                value = varargin{np};
                set(xlh(i), property, value);
            end
        end
        varargout{1} = xlh;
        
    case 'bb'
        if ~isempty(varargin) && all(size(varargin{1})==[2 3]), st{curfig}.bb = varargin{1}; end
        bbox;
        redraw_all;
        
    case 'redraw'
        if nargin < 2
            curfig=GetCurFig;
        else
            curfig=varargin{1};
        end
        redraw_all(curfig);
        eval(st{curfig}.callback);
        if isfield(st{curfig},'registry')
            spm_XYZreg('SetCoords',st{curfig}.centre,st{curfig}.registry.hReg,st{curfig}.registry.hMe);
        end
        
    case 'reload_mats'
        if nargin > 1
            handles = valid_handles(varargin{1});
        else
            handles = valid_handles;
        end
        for i = handles
            fnm = spm_file(st{curfig}.vols{i}.fname, 'number', st{curfig}.vols{i}.n);
            st{curfig}.vols{i}.mat = spm_get_space(fnm);
        end
        % redraw_all (done in w_spm_orthviews('reorient','context_quit'))
        
    case 'reposition'
        if isempty(varargin), tmp = findcent;
        else tmp = varargin{1}; end
        if nargin > 2
            curfig=varargin{2};
        end
        if numel(tmp) == 3
            h = valid_handles(st{curfig}.snap);
            if ~isempty(h)
                tmp = st{curfig}.vols{h(1)}.mat * ...
                    round(st{curfig}.vols{h(1)}.mat\[tmp(:); 1]);
            end
            st{curfig}.centre = tmp(1:3);
        end
        if st{curfig}.yoke
            allfig=curfig;
            for fig=1:numel(st)
                if fig~=curfig
                    if ~isempty(st{fig}) && st{fig}.yoke
                        allfig=[allfig,fig];
                        st{fig}.centre = tmp(1:3);
                    end
                end
            end
        else
            allfig=curfig;
        end
        for fig=allfig
            redraw_all(fig);
            eval(st{fig}.callback);
            if isfield(st{fig},'registry')
                spm_XYZreg('SetCoords',st{fig}.centre,st{fig}.registry.hReg,st{fig}.registry.hMe);
            end
            cm_pos(fig);
        end
        
    case 'setcoords'
        st{curfig}.centre = varargin{1};
        st{curfig}.centre = st{curfig}.centre(:);
        redraw_all;
        eval(st{curfig}.callback);
        cm_pos;
        
    case 'space'
        if numel(varargin) < 1
            st{curfig}.Space = eye(4);
            st{curfig}.bb = maxbb;
            resolution;
            bbox;
            redraw_all;
        else
            space(varargin{:});
            resolution;
            bbox;
            redraw_all;
        end
        
    case 'maxbb'
        st{curfig}.bb = maxbb;
        bbox;
        redraw_all;
        
    case 'resolution'
        resolution(varargin{:});
        bbox;
        redraw_all;
        
    case 'window'
        if numel(varargin)<2
            win = 'auto';
        elseif numel(varargin{2})==2
            win = varargin{2};
        end
        for i=valid_handles(varargin{1})
            st{curfig}.vols{i}.window = win;
        end
        redraw(varargin{1});
        
    case 'delete'
        my_delete(varargin{1});
        
    case 'move'
        move(varargin{1},varargin{2});
        % redraw_all;
        
    case 'reset'
        my_reset;
        
    case 'pos'
        if isempty(varargin)
            H = st{curfig}.centre(:);
        else
            H = pos(varargin{1});
        end
        varargout{1} = H;
        
    case 'interp'
        st{curfig}.hld = varargin{1};
        redraw_all;
        
    case 'xhairs'
        xhairs(varargin{1});
        
    case 'register'
        register(varargin{1});
        
    case 'addblobs'
        addblobs(varargin{:});
        % redraw(varargin{1});
        
    case 'setblobsmax'
        st{curfig}.vols{varargin{1}}.blobs{varargin{2}}.max = varargin{3};
        w_spm_orthviews('redraw')
        
    case 'addcolouredblobs'
        addcolouredblobs(varargin{:});
        % redraw(varargin{1});
        
    case 'addimage'
        addimage(varargin{1}, varargin{2});
        % redraw(varargin{1});
        
    case 'addcolouredimage'
        addcolouredimage(varargin{1}, varargin{2},varargin{3});
        % redraw(varargin{1});
        
    case 'addtruecolourimage'
        if nargin < 2
            varargin(1) = {1};
        end
        if nargin < 3
            varargin(2) = {spm_select(1, 'image', 'Image with activation signal')};
        end
        if nargin < 4
            actc = [];
            while isempty(actc)
                actc = getcmap(spm_input('Colourmap for activation image', '+1','s'));
            end
            varargin(3) = {actc};
        end
        if nargin < 5
            varargin(4) = {0.4};
        end
        if nargin < 6
            actv = spm_vol(varargin{2});
            varargin(5) = {max([eps maxval(actv)])};
        end
        if nargin < 7
            varargin(6) = {min([0 minval(actv)])};
        end
        
        addtruecolourimage(varargin{1}, varargin{2},varargin{3}, varargin{4}, ...
            varargin{5}, varargin{6});
        % redraw(varargin{1});

    case 'settruecolourimage'
        if nargin < 2
            varargin(1) = {1};
        end
        if nargin < 3
            varargin(2) = {spm_select(1, 'image', 'Image with activation signal')};
        end
        if nargin < 4
            actc = [];
            while isempty(actc)
                actc = getcmap(spm_input('Colourmap for activation image', '+1','s'));
            end
            varargin(3) = {actc};
        end
        if nargin < 5
            varargin(4) = {0.4};
        end
        if nargin < 6
            actv = spm_vol(varargin{2});
            varargin(5) = {max([eps maxval(actv)])};
        end
        if nargin < 7
            varargin(6) = {min([0 minval(actv)])};
        end
        
        if nargin < 8
            varargin(7) = 1;
        end
        
        settruecolourimage(varargin{1}, varargin{2},varargin{3}, varargin{4}, ...
            varargin{5}, varargin{6}, varargin{7});
        
    case 'addcolourbar'
        addcolourbar(varargin{1}, varargin{2});
        
    case 'redrawcolourbar'
        curfig=varargin{1};
        bset=varargin{2};
        
        if ~bset
            colormap(gray(64));
            return;
        end
        
        mn=st{curfig}.vols{1}.blobs{bset}.min;
        mx=st{curfig}.vols{1}.blobs{bset}.max;
        
        csz   = size(st{curfig}.vols{1}.blobs{bset}.colour.cmap);
        cdata = reshape(st{curfig}.vols{1}.blobs{bset}.colour.cmap, [csz(1) 1 csz(2)]);
        
        redraw_colourbar(1, bset, [mn, mx], cdata);
        
    case {'removeblobs','rmblobs'}
        rmblobs(varargin{1});
        % redraw(varargin{1});
        
    case 'addcontext'
        if nargin == 1
            handles = 1:max_img;
        else
            handles = varargin{1};
        end
        addcontexts(handles);
        
    case {'removecontext','rmcontext'}
        if nargin == 1
            handles = 1:max_img;
        else
            handles = varargin{1};
        end
        rmcontexts(handles);
        
    case 'context_menu'
        c_menu(varargin{:});
        
    case 'valid_handles'
        if nargin == 1
            handles = 1:max_img;
        else
            handles = varargin{1};
        end
        varargout{1} = valid_handles(handles);

    case 'zoom'
        zoom_op(varargin{:});
        
    case 'zoommenu'
        if isempty(zoomlist)
            zoomlist = [NaN 0 5    10  20 40 80 Inf];
            reslist  = [1   1 .125 .25 .5 .5 1  1  ];
        end
        if nargin >= 3
            if all(cellfun(@isnumeric,varargin(1:2))) && ...
                    numel(varargin{1})==numel(varargin{2})
                zoomlist = varargin{1}(:);
                reslist  = varargin{2}(:);
            else
                warning('w_spm_orthviews:zoom',...
                        'Invalid zoom or resolution list{curfig}.')
            end
        end
        if nargout > 0
            varargout{1} = zoomlist;
        end
        if nargout > 1
            varargout{2} = reslist;
        end
        
    otherwise
        addonaction = strcmpi(st{curfig}.plugins,action);
        if any(addonaction)
            feval(['w_spm_ov_' st{curfig}.plugins{addonaction}],varargin{:});
        end
end

spm('Pointer','Arrow');


%==========================================================================
% function H = specify_image(img)
%==========================================================================
function H = specify_image(img)
global st
curfig=GetCurFig;

H = [];
if isstruct(img)
    V = img(1);
else
    try
        V = spm_vol(img);
    catch
        fprintf('Can not use image "%s"\n', img);
        return;
    end
end
if numel(V)>1, V=V(1); end

ii = 1;
while ~isempty(st{curfig}.vols{ii}), ii = ii + 1; end

if isfield(st{curfig}.vols{1}, 'ax')
    VField=fieldnames(V);
    for n=1:size(VField, 1)
        st{curfig}.vols{1}.(VField{n})=V.(VField{n});
    end
    st{curfig}.bb=[];
    H=1;
    return;
end

FullMatlabVersion = sscanf(version,'%d.%d.%d.%d%s');
DeleteFcn=[];

V.ax = cell(3,1);
handles=guidata(curfig);
for i=1:3
    ax = axes('Visible','on', 'Parent',st{curfig}.fig, ... % 'DrawMode','fast', 
        'YDir','normal', 'DeleteFcn',DeleteFcn, 'ButtonDownFcn',@repos_start,...
        'Parent', handles.ViewFrame);
    d  = image(0, 'Tag','Transverse', 'Parent',ax, 'DeleteFcn',DeleteFcn);
    set(ax, 'Ydir','normal', 'ButtonDownFcn',@repos_start);
    
    lx = line(0,0, 'Parent',ax, 'DeleteFcn',DeleteFcn);
    ly = line(0,0, 'Parent',ax, 'DeleteFcn',DeleteFcn);
    if ~st{curfig}.xhairs
        set(lx, 'Visible','off');
        set(ly, 'Visible','off');
    end
    axis(ax, 'image');
    V.ax{i} = struct('ax',ax,'d',d,'lx',lx,'ly',ly);
end
V.premul    = eye(4);
V.window    = 'auto';
V.mapping   = 'linear';
st{curfig}.vols{ii} = V;

H = ii;


%==========================================================================
% function addblobs(handle, xyz, t, mat, name)
%==========================================================================
function addblobs(handle, xyz, t, mat, name)
global st
curfig=GetCurFig;
if nargin < 5
    name = '';
end
for i=valid_handles(handle)
    if ~isempty(xyz)
        rcp      = round(xyz);
        dim      = max(rcp,[],2)';
        off      = rcp(1,:) + dim(1)*(rcp(2,:)-1 + dim(2)*(rcp(3,:)-1));
        vol      = zeros(dim)+NaN;
        vol(off) = t;
        vol      = reshape(vol,dim);
        st{curfig}.vols{i}.blobs=cell(1,1);
        mx = max([eps max(t)]);
        mn = min([0 min(t)]);
        st{curfig}.vols{i}.blobs{1} = struct('vol',vol,'mat',mat,'max',mx, 'min',mn,'name',name);
        addcolourbar(handle,1);
    end
end


%==========================================================================
% function addimage(handle, fname)
%==========================================================================
function addimage(handle, fname)
global st
curfig=GetCurFig;
for i=valid_handles(handle)
    if isstruct(fname)
        vol = fname(1);
    else
        vol = spm_vol(fname);
    end
    mat = vol.mat;
    st{curfig}.vols{i}.blobs=cell(1,1);
    mx = max([eps maxval(vol)]);
    mn = min([0 minval(vol)]);
    st{curfig}.vols{i}.blobs{1} = struct('vol',vol,'mat',mat,'max',mx,'min',mn);
    addcolourbar(handle,1);
end


%==========================================================================
% function addcolouredblobs(handle, xyz, t, mat, colour, name)
%==========================================================================
function addcolouredblobs(handle, xyz, t, mat, colour, name)
if nargin < 6
    name = '';
end
global st
curfig=GetCurFig;
for i=valid_handles(handle)
    if ~isempty(xyz)
        rcp      = round(xyz);
        dim      = max(rcp,[],2)';
        off      = rcp(1,:) + dim(1)*(rcp(2,:)-1 + dim(2)*(rcp(3,:)-1));
        vol      = zeros(dim)+NaN;
        vol(off) = t;
        vol      = reshape(vol,dim);
        if ~isfield(st{curfig}.vols{i},'blobs')
            st{curfig}.vols{i}.blobs=cell(1,1);
            bset = 1;
        else
            bset = numel(st{curfig}.vols{1}.blobs)+1;
        end
        mx = max([eps maxval(vol)]);
        mn = min([0 minval(vol)]);
        st{curfig}.vols{1}.blobs{bset} = struct('vol',vol, 'mat',mat, ...
            'max',mx, 'min',mn, 'colour',colour, 'name',name);
    end
end


%==========================================================================
% function addcolouredimage(handle, fname,colour)
%==========================================================================
function addcolouredimage(handle, fname,colour)
global st
curfig=GetCurFig;
for i=valid_handles(handle)
    if isstruct(fname)
        vol = fname(1);
    else
        vol = spm_vol(fname);
    end
    mat = vol.mat;
    if ~isfield(st{curfig}.vols{i},'blobs')
        st{curfig}.vols{i}.blobs=cell(1,1);
        bset = 1;
    else
        bset = numel(st{curfig}.vols{i}.blobs)+1;
    end
    mx = max([eps maxval(vol)]);
    mn = min([0 minval(vol)]);
    st{curfig}.vols{i}.blobs{bset} = struct('vol',vol, 'mat',mat, ...
        'max',mx, 'min',mn, 'colour',colour);
end


%==========================================================================
% function addtruecolourimage(handle,fname,colourmap,prop,mx,mn)
%==========================================================================
function addtruecolourimage(curfig,fname,colourmap,prop,mx,mn)
% adds true colour image to current displayed image
global st
%Remove loop by Sandy
%for i=valid_handles(handle)
if isstruct(fname)
    vol = fname(1);
else
    vol = spm_vol(fname);
end
mat = vol.mat;
if ~isfield(st{curfig}.vols{1},'blobs')
    st{curfig}.vols{1}.blobs=cell(1,1);
    bset = 1;
else
    bset = numel(st{curfig}.vols{1}.blobs)+1;
end

c = struct('cmap', colourmap,'prop',prop);
st{curfig}.vols{1}.blobs{bset} = struct('vol',vol, 'mat',mat, ...
    'max',mx, 'min',mn, 'colour',c);
addcolourbar(1,bset);
st{curfig}.curblob=bset;

%Add by Sandy to change colourbar when add overlay
if ~isfield(st{curfig}.vols{1}.blobs{bset},'colour')
    cmap = get(st{curfig}.fig,'Colormap');
    if size(cmap,1)~=128
        figure(st{curfig}.fig)
        spm_figure('Colormap','gray-hot')
    end
    redraw_colourbar(1,bset,[mn mx],(1:64)'+64);
elseif isstruct(st{curfig}.vols{1}.blobs{bset}.colour)
    csz   = size(st{curfig}.vols{1}.blobs{bset}.colour.cmap);
    cdata = reshape(st{curfig}.vols{1}.blobs{bset}.colour.cmap, [csz(1) 1 csz(2)]);
    redraw_colourbar(1,bset,[mn mx],cdata);
end
%end

%==========================================================================
% function settruecolourimage(handle,fname,colourmap,prop,mx,mn)
%==========================================================================
function settruecolourimage(curfig,fname,colourmap,prop,mx,mn,bset)
% set true colour image to current displayed image Add by Sandy
global st

if isstruct(fname)
    vol = fname(1);
else
    vol = spm_vol(fname);
end
mat = vol.mat;

c = struct('cmap', colourmap,'prop',prop);
st{curfig}.vols{1}.blobs{bset} = struct('vol',vol, 'mat',mat, ...
    'max',mx, 'min',mn, 'colour',c);
addcolourbar(1, bset, curfig);
st{curfig}.curblob=bset;
if ~isfield(st{curfig}.vols{1}.blobs{bset},'colour')
    cmap = get(st{curfig}.fig,'Colormap');
    if size(cmap,1)~=128
        figure(st{curfig}.fig)
        spm_figure('Colormap','gray-hot')
    end
    redraw_colourbar(1,bset,[mn mx],(1:64)'+64,curfig);
elseif isstruct(st{curfig}.vols{1}.blobs{bset}.colour)
    csz   = size(st{curfig}.vols{1}.blobs{bset}.colour.cmap);
    cdata = reshape(st{curfig}.vols{1}.blobs{bset}.colour.cmap, [csz(1) 1 csz(2)]);
    redraw_colourbar(1,bset,[mn mx],cdata,curfig);
end
st{curfig}.weight=zeros(size(st{curfig}.vols{1}.blobs{1}.vol.Data));

%==========================================================================
% function addcolourbar(vh,bh)
%==========================================================================
function addcolourbar(vh,bh,curfig)
global st
if nargin < 3 
    curfig=GetCurFig;
end
if st{curfig}.mode == 0,
    axpos = get(st{curfig}.vols{vh}.ax{2}.ax,'Position');
else
    axpos = get(st{curfig}.vols{vh}.ax{1}.ax,'Position');
end
handles=guidata(st{curfig}.fig);
if ~isfield(handles, 'DPABI_fig');
    st{curfig}.vols{vh}.blobs{bh}.cbar = axes('Parent',st{curfig}.fig,...
        'Position',[(axpos(1)+axpos(3)+0.05+(bh-1)*.1) (axpos(2)+0.005) 0.05 (axpos(4)-0.01)],...
        'Box','on', 'YDir','normal', 'XTickLabel',[], 'XTick',[]);
else
    st{curfig}.vols{vh}.blobs{bh}.cbar = handles.ColorAxe;
end
if isfield(st{curfig}.vols{vh}.blobs{bh},'name')
    ylabel(st{curfig}.vols{vh}.blobs{bh}.name,'parent',st{curfig}.vols{vh}.blobs{bh}.cbar);
end

%==========================================================================
% function rmblobs(handle)
%==========================================================================
function rmblobs(handle)
global st
curfig=GetCurFig;
%Remove loop by Sandy
if ~isempty(st{curfig})
    if isfield(st{curfig}.vols{1},'blobs')
        handles=guidata(st{curfig}.fig);
        if ~isfield(handles, 'DPABI_fig');
            for j=1:numel(st{curfig}.vols{1}.blobs)
                if isfield(st{curfig}.vols{1}.blobs{j},'cbar') && ishandle(st{curfig}.vols{1}.blobs{j}.cbar),
                    delete(st{curfig}.vols{1}.blobs{j}.cbar);
                end
            end
        end
        st{curfig}.vols{1} = rmfield(st{curfig}.vols{1},'blobs');
    end
    st{curfig}.curblob=0;
end


%==========================================================================
% function register(hreg)
%==========================================================================
function register(hreg)
global st
curfig=GetCurFig;
%tmp = uicontrol('Position',[0 0 1 1],'Visible','off','Parent',st{curfig}.fig);
h   = valid_handles;
if ~isempty(h)
    tmp = st{curfig}.vols{h(1)}.ax{1}.ax;
    st{curfig}.registry = struct('hReg',hreg,'hMe', tmp);
    spm_XYZreg('Add2Reg',st{curfig}.registry.hReg,st{curfig}.registry.hMe, 'w_spm_orthviews');
else
    warning('Nothing to register with');
end
st{curfig}.centre = spm_XYZreg('GetCoords',st{curfig}.registry.hReg);
st{curfig}.centre = st{curfig}.centre(:);


%==========================================================================
% function xhairs(state)
%==========================================================================
function xhairs(state)
global st
curfig=GetCurFig;
st{curfig}.xhairs = 0;
opt = 'on';
if ~strcmpi(state,'on')
    opt = 'off';
else
    st{curfig}.xhairs = 1;
end
for i=valid_handles
    for j=1:3
        set(st{curfig}.vols{i}.ax{j}.lx,'Visible',opt);
        set(st{curfig}.vols{i}.ax{j}.ly,'Visible',opt);
    end
end


%==========================================================================
% function H = pos(handle)
%==========================================================================
function H = pos(handle)
global st
curfig=GetCurFig;
H = [];
for i=valid_handles(handle)
    is = inv(st{curfig}.vols{i}.premul*st{curfig}.vols{i}.mat);
    H = is(1:3,1:3)*st{curfig}.centre(:) + is(1:3,4);
end


%==========================================================================
% function my_reset
%==========================================================================
function my_reset
global st
curfig=GetCurFig;
if ~isempty(st{curfig}) && isfield(st{curfig},'registry') && ishandle(st{curfig}.registry.hMe)
    delete(st{curfig}.registry.hMe); st = rmfield(st,'registry');
end
my_delete(curfig);
reset_st;


%==========================================================================
% function my_delete(handle)
%==========================================================================
function my_delete(handle)
global st
curfig=GetCurFig;
% remove blobs (and colourbars, if any)
rmblobs(handle);
% remove displayed axes
% Remove loop by Sandy
MainHandle=guidata(curfig);
kids = get(MainHandle.ViewFrame,'Children');
for j=1:3
    try
        if any(kids == st{curfig}.vols{1}.ax{j}.ax)
            set(get(st{curfig}.vols{1}.ax{j}.ax,'Children'),'DeleteFcn','');
            delete(st{curfig}.vols{1}.ax{j}.ax);
        end
    end
end
st{curfig}.vols{1} = [];
st{curfig}=[];


%==========================================================================
% function resolution(res)
%==========================================================================
function resolution(res)
global st
curfig=GetCurFig;
if ~nargin, res = 1; end % Default minimum resolution 1mm
for i=valid_handles
    % adapt resolution to smallest voxel size of displayed images
    res  = min([res,sqrt(sum((st{curfig}.vols{i}.mat(1:3,1:3)).^2))]);
end
res      = res/mean(svd(st{curfig}.Space(1:3,1:3)));
Mat      = diag([res res res 1]);
st{curfig}.Space = st{curfig}.Space*Mat;
st{curfig}.bb    = st{curfig}.bb/res;


%==========================================================================
% function move(handle,pos)
%==========================================================================
function move(handle,pos)
global st
curfig=GetCurFig;
for i=valid_handles(handle)
    st{curfig}.vols{i}.area = pos;
end
bbox;
% redraw(valid_handles(handle));


%==========================================================================
% function bb = maxbb
%==========================================================================
function bb = maxbb
global st
curfig=GetCurFig;
mn = [Inf Inf Inf];
mx = -mn;
for i=valid_handles
    premul = st{curfig}.Space \ st{curfig}.vols{i}.premul;
    bb = spm_get_bbox(st{curfig}.vols{i}, 'fv', premul);
    mx = max([bb ; mx]);
    mn = min([bb ; mn]);
end
bb = [mn ; mx];


%==========================================================================
% function space(handle,M,dim)
%==========================================================================
function space(handle,M,dim)
global st
curfig=GetCurFig;
if ~isempty(st{curfig}.vols{handle})
    if nargin < 2
        M = st{curfig}.vols{handle}.mat;
        dim = st{curfig}.vols{handle}.dim(1:3);
    end
    Mat   = st{curfig}.vols{handle}.premul(1:3,1:3)*M(1:3,1:3);
    vox   = sqrt(sum(Mat.^2));
    if det(Mat(1:3,1:3))<0, vox(1) = -vox(1); end
    Mat   = diag([vox 1]);
    Space = (M)/Mat;
    bb    = [1 1 1; dim];
    bb    = [bb [1;1]];
    bb    = bb*Mat';
    bb    = bb(:,1:3);
    bb    = sort(bb);
    st{curfig}.Space = Space;
    st{curfig}.bb = bb;
end


%==========================================================================
% function zoom_op(fov,res)
%==========================================================================
function zoom_op(fov,res)
global st
curfig=GetCurFig;
if nargin < 1, fov = Inf; end
if nargin < 2, res = Inf; end

if isinf(fov)
    st{curfig}.bb = maxbb;
elseif isnan(fov) || fov == 0
    current_handle = valid_handles;
    if numel(current_handle) > 1 % called from check reg context menu
        current_handle = get_current_handle;
    end
    if fov == 0
        % zoom to bounding box of current image ~= 0
        thr = 'nz';
    else
        % zoom to bounding box of current image > chosen threshold
        thr = spm_input('Threshold (Y > ...)', '+1', 'r', '0', 1);
    end
    premul = st{curfig}.Space \ st{curfig}.vols{current_handle}.premul;
    st{curfig}.bb = spm_get_bbox(st{curfig}.vols{current_handle}, thr, premul);
else
    vx    = sqrt(sum(st{curfig}.Space(1:3,1:3).^2));
    vx    = vx.^(-1);
    pos   = w_spm_orthviews('pos');
    pos   = st{curfig}.Space\[pos ; 1];
    pos   = pos(1:3)';
    st{curfig}.bb = [pos-fov*vx; pos+fov*vx];
end
resolution(res);
bbox;
redraw_all;
if isfield(st{curfig}.vols{1},'sdip')
    spm_eeg_inv_vbecd_disp('RedrawDip');
end


%==========================================================================
% function repos_start(varargin)
% function repos_move(varargin)
% function repos_end(varargin)
%==========================================================================
function repos_start(varargin)
% don't use right mouse button to start reposition
if ~strcmpi(get(gcbf,'SelectionType'),'alt')
    set(gcbf,'windowbuttonmotionfcn',@repos_move, 'windowbuttonupfcn',@repos_end);
    w_spm_orthviews('reposition');
end

function repos_move(varargin)
w_spm_orthviews('reposition');

function repos_end(varargin)
set(gcbf,'windowbuttonmotionfcn','', 'windowbuttonupfcn','');


%==========================================================================
% function bbox
%==========================================================================
function bbox
global st
curfig=GetCurFig;
Dims = diff(st{curfig}.bb)'+1;

TD = Dims([1 2])';
CD = Dims([1 3])';
if st{curfig}.mode == 0, SD = Dims([3 2])'; else SD = Dims([2 3])'; end
%Add by Sandy for DPABI_VIEW call
handles=guidata(st{curfig}.fig);
if ~isfield(handles, 'DPABI_fig');
    un    = get(st{curfig}.fig,'Units');set(st{curfig}.fig,'Units','Pixels');
    sz    = get(st{curfig}.fig,'Position');set(st{curfig}.fig,'Units',un);
    sz    = sz(3:4);
    sz(1) = sz(1);
    sz(2) = sz(2);
else
    un    = get(handles.ViewFrame,'Units');set(handles.ViewFrame,'Units','Pixels');
    sz    = get(handles.ViewFrame,'Position');set(handles.ViewFrame,'Units',un);
    
    offxy = sz(1:2);
    sz    = sz(3:4);
end

for i=valid_handles
    area   = st{curfig}.vols{i}.area(:);
    area   = [area(1)*sz(1) area(2)*sz(2) area(3)*sz(1) area(4)*sz(2)];
    if st{curfig}.mode == 0
        sx = area(3)/(Dims(1)+Dims(3))/1.02;
    else
        sx = area(3)/(Dims(1)+Dims(2))/1.02;
    end
    sy     = area(4)/(Dims(2)+Dims(3))/1.02;
    s      = min([sx sy]);
    
    offy   = (area(4)-(Dims(2)+Dims(3))*1.02*s)/2 + area(2);
    sky    = s*(Dims(2)+Dims(3))*0.02;
    if st{curfig}.mode == 0
        offx = (area(3)-(Dims(1)+Dims(3))*1.02*s)/2 + area(1);
        skx  = s*(Dims(1)+Dims(3))*0.02;
    else
        offx = (area(3)-(Dims(1)+Dims(2))*1.02*s)/2 + area(1);
        skx  = s*(Dims(1)+Dims(2))*0.02;
    end
    
    %if isfield(handles, 'DPABI_fig');
    %    offx=offx+offxy(1);
    %    offy=offy+offxy(2);
    %end
    
    % Transverse
    set(st{curfig}.vols{i}.ax{1}.ax,'Units','pixels', ...
        'Position',[offx offy s*Dims(1) s*Dims(2)],...
        'Units','normalized','Xlim',[0 TD(1)]+0.5,'Ylim',[0 TD(2)]+0.5,...
        'Visible','on','XTick',[],'YTick',[]);
    
    % Coronal
    set(st{curfig}.vols{i}.ax{2}.ax,'Units','Pixels',...
        'Position',[offx offy+s*Dims(2)+sky s*Dims(1) s*Dims(3)],...
        'Units','normalized','Xlim',[0 CD(1)]+0.5,'Ylim',[0 CD(2)]+0.5,...
        'Visible','on','XTick',[],'YTick',[]);
    
    % Sagittal
    if st{curfig}.mode == 0
        set(st{curfig}.vols{i}.ax{3}.ax,'Units','Pixels', 'Box','on',...
            'Position',[offx+s*Dims(1)+skx offy s*Dims(3) s*Dims(2)],...
            'Units','normalized','Xlim',[0 SD(1)]+0.5,'Ylim',[0 SD(2)]+0.5,...
            'Visible','on','XTick',[],'YTick',[]);
    else
        set(st{curfig}.vols{i}.ax{3}.ax,'Units','Pixels', 'Box','on',...
            'Position',[offx+s*Dims(1)+skx offy+s*Dims(2)+sky s*Dims(2) s*Dims(3)],...
            'Units','normalized','Xlim',[0 SD(1)]+0.5,'Ylim',[0 SD(2)]+0.5,...
            'Visible','on','XTick',[],'YTick',[]);
    end
end


%==========================================================================
% function mx = maxval(vol)
%==========================================================================
function mx = maxval(vol)
if isstruct(vol)
    mx = -Inf;
    for i=1:vol.dim(3)
        
        if ~isfield(vol,'Data')
            tmp = spm_slice_vol(vol,spm_matrix([0 0 i]),vol.dim(1:2),0);
        else   %Revised by YAN Chao-Gan, 130720. Could also work with Data has been read into memory other than only depending on the file.
            tmp = spm_slice_vol(vol.Data,spm_matrix([0 0 i]),vol.dim(1:2),0);
        end
        
        
        imx = max(tmp(isfinite(tmp)));
        if ~isempty(imx), mx = max(mx,imx); end
    end
else
    mx = max(vol(isfinite(vol)));
end


%==========================================================================
% function mn = minval(vol)
%==========================================================================
function mn = minval(vol)
if isstruct(vol)
    mn = Inf;
    for i=1:vol.dim(3)
        
        if ~isfield(vol,'Data')
            tmp = spm_slice_vol(vol,spm_matrix([0 0 i]),vol.dim(1:2),0);
        else   %Revised by YAN Chao-Gan, 130720. Could also work with Data has been read into memory other than only depending on the file.
            tmp = spm_slice_vol(vol.Data,spm_matrix([0 0 i]),vol.dim(1:2),0);
        end
        
        
        imn = min(tmp(isfinite(tmp)));
        if ~isempty(imn), mn = min(mn,imn); end
    end
else
    mn = min(vol(isfinite(vol)));
end


%==========================================================================
% function redraw(arg1)
%==========================================================================
function redraw(fig)
global st
for curfig=fig

bb   = st{curfig}.bb;
%Add by Sandy for DPABI_VIEW call

handles=guidata(st{curfig}.fig);
if isfield(handles, 'DPABI_fig')
    X=st{curfig}.centre(1);
    if X<=bb(1,1)
        X=bb(1,1);
        set(handles.XReduceBtn, 'Enable', 'Off');
        set(handles.IAddBtn, 'Enable', 'Off');
    else
        set(handles.XReduceBtn, 'Enable', 'On');
        set(handles.IAddBtn, 'Enable', 'On');
    end
    if X>=bb(2,1)
        X=bb(2,1);
        set(handles.XAddBtn, 'Enable', 'Off');
        set(handles.IReduceBtn, 'Enable', 'Off');
    else
        set(handles.XAddBtn, 'Enable', 'On');
        set(handles.IReduceBtn, 'Enable', 'On');
    end
    
    Y=st{curfig}.centre(2);
    if Y<=bb(1,2)
        Y=bb(1,2);
        set(handles.YReduceBtn, 'Enable', 'Off');
        set(handles.JReduceBtn, 'Enable', 'Off');
    else
        set(handles.YReduceBtn, 'Enable', 'On');
        set(handles.JReduceBtn, 'Enable', 'On');
    end
    if Y>=bb(2,2)
        Y=bb(2,2);
        set(handles.YAddBtn, 'Enable', 'Off');
        set(handles.JAddBtn, 'Enable', 'Off');
    else
        set(handles.YAddBtn, 'Enable', 'On');
        set(handles.JAddBtn, 'Enable', 'On');
    end
    
    Z=st{curfig}.centre(3); 
    if Z<=bb(1,3)
        Z=bb(1,3);
        set(handles.ZReduceBtn, 'Enable', 'Off');
        set(handles.KReduceBtn, 'Enable', 'Off');
    else
        set(handles.ZReduceBtn, 'Enable', 'On');
        set(handles.KReduceBtn, 'Enable', 'On');
    end
    if Z>=bb(2,3)
        Z=bb(2,3);
        set(handles.ZAddBtn, 'Enable', 'Off');
        set(handles.KAddBtn, 'Enable', 'Off');
    else
        set(handles.ZAddBtn, 'Enable', 'On');
        set(handles.KAddBtn, 'Enable', 'On');
    end
    
    set(handles.XEntry,...
        'String', sprintf('%d', round(X)));
    set(handles.YEntry,...
        'String', sprintf('%d', round(Y)));
    set(handles.ZEntry,...
        'String', sprintf('%d', round(Z)));
    
    tmp=inv(st{curfig}.vols{1}.mat)*[X;Y;Z;1];
    I=round(tmp(1));
    J=round(tmp(2));
    K=round(tmp(3));
    
    AString={};
    if ~isempty(st{curfig}.AtlasInfo)
        for idx=1:numel(st{curfig}.AtlasInfo)
            AStruct=st{curfig}.AtlasInfo{idx};
            AMat=AStruct.Template.mat;
            APos=round(inv(AMat)*[X;Y;Z;1]);
            AI=APos(1);
            AJ=APos(2);
            AK=APos(3);
            try
                AIndex=AStruct.Template.Data(AI, AJ, AK);
            catch
                AIndex=0;
            end
            AName=AStruct.Template.Alias;
            ALab =cellfun(@(x) isequal(x, AIndex),...
                AStruct.Reference(:, 2));
            if any(ALab)
                ARegion=AStruct.Reference{ALab, 1};
            else
                ARegion='None';
            end
            AString{numel(AString)+1, 1}=...
                sprintf(' (%s) %s', AName, ARegion);
        end
    end
    
    if isfield(st{curfig}.vols{1}, 'blobs')
        curblob=st{curfig}.curblob;
        tmp=inv(st{curfig}.vols{1}.blobs{curblob}.vol.mat)*[X;Y;Z;1];
        OI=round(tmp(1));
        OJ=round(tmp(2));
        OK=round(tmp(3));

        set(handles.IEntry, 'String', sprintf('%d', OI));
        set(handles.JEntry, 'String', sprintf('%d', OJ));
        set(handles.KEntry, 'String', sprintf('%d', OK));
        try 
            OverlayValue=st{curfig}.vols{1}.blobs{curblob}.vol.Data(OI,OJ,OK);
        catch
            OverlayValue=0;
        end
        set(handles.OverlayValue, 'String',...
            sprintf('%g', OverlayValue));
        TCFlag=st{curfig}.TCFlag;
        if ~isempty(TCFlag)
            if ~ishandle(TCFlag)
                st{curfig}.TCFlag=0;
                if isfield(st{curfig}, 'TCLinObj')
                    st{curfig}=rmfield(st{curfig}, 'TCLinObj');
                end
            else
                TCHandle=guidata(TCFlag);
                Headers=TCHandle.Headers;
                TCAxe=TCHandle.TCAxe;
                if ~isfield(st{curfig}, 'TCLinObj') || ...
                    isempty(st{curfig}.TCLinObj{1}) || ...    
                    ~ishandle(st{curfig}.TCLinObj{1})
                    st{curfig}.TCLinObj=cell(numel(Headers), 1);
                end
                set(TCHandle.CoordFrame, 'Title',...
                    sprintf('Coordinate: X=%+d, Y=%+d, Z=%+d; I=%d, J=%d, K=%d',...
                    round(X), round(Y), round(Z), OI, OJ, OK))
                
                for i=1:numel(Headers)
                    if ~isempty(Headers{i})
                        tmp=inv(Headers{i}.mat)*[X;Y;Z;1];
                        TCI=round(tmp(1));
                        TCJ=round(tmp(2));
                        TCK=round(tmp(3));
                        try
                            TC=squeeze(Headers{i}.Raw(TCI,TCJ,TCK,:));
                        catch
                            TC=zeros(size(Headers{i}.Raw, 4), 1);
                        end
                        
                        if get(TCHandle.FrequencyButton, 'Value')
                            TR=str2double(get(TCHandle.TREntry, 'String'));
                            TNum=length(TC);
                            Pad=2^nextpow2(TNum);
                            
                            TC=TC-mean(TC);
                            AM=2*abs(fft(TC, Pad))/TNum;
                            AMP=(1:length(AM))';
                            AMP=(AMP-1)/(Pad*TR);
                            TC=AM(1:ceil(length(AM)/2)+1);
                            TCP=AMP(1:ceil(length(AMP)/2)+1);
                        else
                            TCP=(1:size(Headers{i}.Raw, 4))';
                        end
                        
                        if ~isempty(st{curfig}.TCLinObj{i})
                            set(st{curfig}.TCLinObj{i}, 'XData', TCP);
                            set(st{curfig}.TCLinObj{i}, 'YData', TC);
                        else
                            lin_obj=plot(TCAxe, TCP, TC);
                            legend(lin_obj, Headers{i}.fname, ...
                                'Location', 'NorthOutside');
                            st{curfig}.TCLinObj{i}=lin_obj;
                        end
                    end
                end
            end
        end
    else
        set(handles.IEntry, 'String', sprintf('%d', I));
        set(handles.JEntry, 'String', sprintf('%d', J));
        set(handles.KEntry, 'String', sprintf('%d', K));        
        set(handles.OverlayValue, 'String', '0');
        OverlayValue=0;
    end
    try
        UnderlayValue=st{curfig}.vols{1}.Data(I,J,K);
    catch
        UnderlayValue=0;
    end
    set(handles.UnderlayValue, 'String',...
        sprintf('%g', UnderlayValue));
    
    UOString={'';...
        sprintf(' Underlay Value -> %g', UnderlayValue);...
        sprintf(' Overlay Value -> %g', OverlayValue);...
        ''};
    set(handles.AtlasEntry, 'String', [UOString; AString]);
    
    MPFlag=st{curfig}.MPFlag;
    if ~isempty(MPFlag)
        for m=1:numel(MPFlag)
            if ~isempty(st{curfig}.MPFlag{m}) && ishandle(st{curfig}.MPFlag{m})
                w_Montage('RedrawXhairs', curfig, MPFlag{m});
            end
        end
    end
end

Dims = round(diff(bb)'+1);
is   = inv(st{curfig}.Space);
cent = is(1:3,1:3)*st{curfig}.centre(:) + is(1:3,4);

%for i = valid_handles(arg1) Remove loop by Sandy
    i=1;
    M = st{curfig}.Space\st{curfig}.vols{i}.premul*st{curfig}.vols{i}.mat;
    TM0 = [ 1 0 0 -bb(1,1)+1
            0 1 0 -bb(1,2)+1
            0 0 1 -cent(3)
            0 0 0 1];
    TM = inv(TM0*M);
    TD = Dims([1 2]);
    
    CM0 = [ 1 0 0 -bb(1,1)+1
            0 0 1 -bb(1,3)+1
            0 1 0 -cent(2)
            0 0 0 1];
    CM = inv(CM0*M);
    CD = Dims([1 3]);
    
    if st{curfig}.mode ==0
        SM0 = [ 0 0 1 -bb(1,3)+1
                0 1 0 -bb(1,2)+1
                1 0 0 -cent(1)
                0 0 0 1];
        SM = inv(SM0*M); 
        SD = Dims([3 2]);
    else
        SM0 = [ 0 -1 0 +bb(2,2)+1
                0  0 1 -bb(1,3)+1
                1  0 0 -cent(1)
                0  0 0 1];
        SM = inv(SM0*M);
        SD = Dims([2 3]);
    end
    
    try
        
        %Revised by YAN Chao-Gan, 130720. Could also work with Data has been read into memory other than only depending on the file.
        if ~isfield(st{curfig}.vols{i},'Data')
            imgt = spm_slice_vol(st{curfig}.vols{i},TM,TD,st{curfig}.hld)';
            imgc = spm_slice_vol(st{curfig}.vols{i},CM,CD,st{curfig}.hld)';
            imgs = spm_slice_vol(st{curfig}.vols{i},SM,SD,st{curfig}.hld)';
        else   %Revised by YAN Chao-Gan, 130720. Could also work with Data has been read into memory other than only depending on the file.
            imgt = spm_slice_vol(st{curfig}.vols{i}.Data,TM,TD,st{curfig}.hld)';
            imgc = spm_slice_vol(st{curfig}.vols{i}.Data,CM,CD,st{curfig}.hld)';
            imgs = spm_slice_vol(st{curfig}.vols{i}.Data,SM,SD,st{curfig}.hld)';
        end
        
        ok   = true;
    catch
        fprintf('Cannot access file "%s".\n', st{curfig}.vols{i}.fname);
        fprintf('%s\n',getfield(lasterror,'message'));
        ok   = false;
    end
    if ok
        % get min/max threshold
        if strcmp(st{curfig}.vols{i}.window,'auto')
            mn = -Inf;
            mx = Inf;
        else
            mn = min(st{curfig}.vols{i}.window);
            mx = max(st{curfig}.vols{i}.window);
        end
        % threshold images
        imgt = max(imgt,mn); imgt = min(imgt,mx);
        imgc = max(imgc,mn); imgc = min(imgc,mx);
        imgs = max(imgs,mn); imgs = min(imgs,mx);
        % compute intensity mapping, if histeq is available
        if license('test','image_toolbox') == 0
            st{curfig}.vols{i}.mapping = 'linear';
        end
        switch st{curfig}.vols{i}.mapping
            case 'linear'
            case 'histeq'
                % scale images to a range between 0 and 1
                imgt1=(imgt-min(imgt(:)))/(max(imgt(:)-min(imgt(:)))+eps);
                imgc1=(imgc-min(imgc(:)))/(max(imgc(:)-min(imgc(:)))+eps);
                imgs1=(imgs-min(imgs(:)))/(max(imgs(:)-min(imgs(:)))+eps);
                img  = histeq([imgt1(:); imgc1(:); imgs1(:)],1024);
                imgt = reshape(img(1:numel(imgt1)),size(imgt1));
                imgc = reshape(img(numel(imgt1)+(1:numel(imgc1))),size(imgc1));
                imgs = reshape(img(numel(imgt1)+numel(imgc1)+(1:numel(imgs1))),size(imgs1));
                mn = 0;
                mx = 1;
            case 'quadhisteq'
                % scale images to a range between 0 and 1
                imgt1=(imgt-min(imgt(:)))/(max(imgt(:)-min(imgt(:)))+eps);
                imgc1=(imgc-min(imgc(:)))/(max(imgc(:)-min(imgc(:)))+eps);
                imgs1=(imgs-min(imgs(:)))/(max(imgs(:)-min(imgs(:)))+eps);
                img  = histeq([imgt1(:).^2; imgc1(:).^2; imgs1(:).^2],1024);
                imgt = reshape(img(1:numel(imgt1)),size(imgt1));
                imgc = reshape(img(numel(imgt1)+(1:numel(imgc1))),size(imgc1));
                imgs = reshape(img(numel(imgt1)+numel(imgc1)+(1:numel(imgs1))),size(imgs1));
                mn = 0;
                mx = 1;
            case 'loghisteq'
                sw = warning('off','MATLAB:log:logOfZero');
                imgt = log(imgt-min(imgt(:)));
                imgc = log(imgc-min(imgc(:)));
                imgs = log(imgs-min(imgs(:)));
                warning(sw);
                imgt(~isfinite(imgt)) = 0;
                imgc(~isfinite(imgc)) = 0;
                imgs(~isfinite(imgs)) = 0;
                % scale log images to a range between 0 and 1
                imgt1=(imgt-min(imgt(:)))/(max(imgt(:)-min(imgt(:)))+eps);
                imgc1=(imgc-min(imgc(:)))/(max(imgc(:)-min(imgc(:)))+eps);
                imgs1=(imgs-min(imgs(:)))/(max(imgs(:)-min(imgs(:)))+eps);
                img  = histeq([imgt1(:); imgc1(:); imgs1(:)],1024);
                imgt = reshape(img(1:numel(imgt1)),size(imgt1));
                imgc = reshape(img(numel(imgt1)+(1:numel(imgc1))),size(imgc1));
                imgs = reshape(img(numel(imgt1)+numel(imgc1)+(1:numel(imgs1))),size(imgs1));
                mn = 0;
                mx = 1;
        end
        % recompute min/max for display
        if strcmp(st{curfig}.vols{i}.window,'auto')
            mx = -inf; mn = inf;
        end
        %Add by Sandy, make the same mn/mx in a volume
        if ~isfield(st{curfig}.vols{i},'Data')
            if ~isempty(imgt)
                tmp = imgt(isfinite(imgt));
                mx = max([mx max(max(tmp))]);
            mn = min([mn min(min(tmp))]);
            end
            if ~isempty(imgc)
                tmp = imgc(isfinite(imgc));
                mx = max([mx max(max(tmp))]);
                mn = min([mn min(min(tmp))]);
            end
            if ~isempty(imgs)
                tmp = imgs(isfinite(imgs));
                mx = max([mx max(max(tmp))]);
                mn = min([mn min(min(tmp))]);
            end
        else
            mx=max(max(max(st{curfig}.vols{i}.Data)));
            mn=min(min(min(st{curfig}.vols{i}.Data)));
        end
        if mx==mn, mx=mn+eps; end
        
        if isfield(st{curfig}.vols{i},'blobs')
            j=1;
            if ~isfield(st{curfig}.vols{i}.blobs{j},'colour')
                % Add blobs for display using the split colourmap
                scal = 64/(mx-mn);
                dcoff = -mn*scal;
                imgt = imgt*scal+dcoff;
                imgc = imgc*scal+dcoff;
                imgs = imgs*scal+dcoff;
                
                if isfield(st{curfig}.vols{i}.blobs{j},'max')
                    mx = st{curfig}.vols{i}.blobs{j}.max;
                else
                    mx = max([eps maxval(st{curfig}.vols{i}.blobs{j}.vol)]);
                    st{curfig}.vols{i}.blobs{j}.max = mx;
                end
                if isfield(st{curfig}.vols{i}.blobs{j},'min')
                    mn = st{curfig}.vols{i}.blobs{j}.min;
                else
                    mn = min([0 minval(st{curfig}.vols{i}.blobs{j}.vol)]);
                    st{curfig}.vols{i}.blobs{j}.min = mn;
                end
                
                vol  = st{curfig}.vols{i}.blobs{j}.vol;
                M    = st{curfig}.Space\st{curfig}.vols{i}.premul*st{curfig}.vols{i}.blobs{j}.mat;
                
                %Revised by YAN Chao-Gan, 130720. Could also work with Data has been read into memory other than only depending on the file.
                if ~isfield(vol,'Data')
                    tmpt = spm_slice_vol(vol,inv(TM0*M),TD,[0 NaN])';
                    tmpc = spm_slice_vol(vol,inv(CM0*M),CD,[0 NaN])';
                    tmps = spm_slice_vol(vol,inv(SM0*M),SD,[0 NaN])';
                else   %Revised by YAN Chao-Gan, 130720. Could also work with Data has been read into memory other than only depending on the file.
                    tmpt = spm_slice_vol(vol.Data,inv(TM0*M),TD,[0 NaN])';
                    tmpc = spm_slice_vol(vol.Data,inv(CM0*M),CD,[0 NaN])';
                    tmps = spm_slice_vol(vol.Data,inv(SM0*M),SD,[0 NaN])';
                end
                
                %tmpt_z = find(tmpt==0);tmpt(tmpt_z) = NaN;
                %tmpc_z = find(tmpc==0);tmpc(tmpc_z) = NaN;
                %tmps_z = find(tmps==0);tmps(tmps_z) = NaN;
                
                sc   = 64/(mx-mn);
                off  = 65.51-mn*sc;
                msk  = find(isfinite(tmpt)); imgt(msk) = off+tmpt(msk)*sc;
                msk  = find(isfinite(tmpc)); imgc(msk) = off+tmpc(msk)*sc;
                msk  = find(isfinite(tmps)); imgs(msk) = off+tmps(msk)*sc;
                
            elseif isstruct(st{curfig}.vols{i}.blobs{j}.colour)
                % Add blobs for display using a defined colourmap
                
                % colourmaps
                gryc = (0:63)'*ones(1,3)/63;
                
                % scale grayscale image, not isfinite -> black
                imgt = scaletocmap(imgt,mn,mx,gryc,65);
                imgc = scaletocmap(imgc,mn,mx,gryc,65);
                imgs = scaletocmap(imgs,mn,mx,gryc,65);
                gryc = [gryc; 0 0 0];

               
                %Modified by Sandy for Multi-Overlay 20140104
                umpt=ones(size(imgt));
                umpc=ones(size(imgc));
                umps=ones(size(imgs));
                umpt=repmat(umpt(:),1,3);
                umpc=repmat(umpc(:),1,3);
                umps=repmat(umps(:),1,3);
                ompt=zeros(size(umpt));
                ompc=zeros(size(umpc));
                omps=zeros(size(umps));
                
                if isfield(handles, 'DPABI_fig')
                    IsPreview=get(handles.PreviewBtn, 'Value');
                    if IsPreview
                        AtlasIdx=get(handles.StructuralPopup, 'Value');
                        AStruct=st{curfig}.AtlasInfo{AtlasIdx};
                        vol=AStruct.Template;
                        M=st{curfig}.Space\st{curfig}.vols{i}.premul*vol.mat;

                        AMat=AStruct.Template.mat;
                        APos=round(inv(AMat)*[X;Y;Z;1]);
                        AI=APos(1);
                        AJ=APos(2);
                        AK=APos(3);
            
                        AIndex=AStruct.Template.Data(AI, AJ, AK);
                                
                        if AIndex
                            SSactc=AStruct.Template.CMap;
                            SSactp=0.8;
                               
                            Mask=double(AStruct.Template.Data==AIndex);
                            SStopc = size(SSactc,1)+1;
                              
                            tmpt = spm_slice_vol(Mask,inv(TM0*M),TD,[0 NaN])';
                            tmpc = spm_slice_vol(Mask,inv(CM0*M),CD,[0 NaN])';
                            tmps = spm_slice_vol(Mask,inv(SM0*M),SD,[0 NaN])';
                   
                            tmpt = scaletocmap(tmpt, 0, 1, SSactc, SStopc);
                            tmpc = scaletocmap(tmpc, 0, 1, SSactc, SStopc);
                            tmps = scaletocmap(tmps, 0, 1, SSactc, SStopc);
                                
                            jmpt = SSactp*repmat((tmpt(:)~=SStopc),1,3);
                            jmpc = SSactp*repmat((tmpc(:)~=SStopc),1,3);
                            jmps = SSactp*repmat((tmps(:)~=SStopc),1,3);
                              
                            umpt = umpt-jmpt;
                            umpc = umpc-jmpc;
                            umps = umps-jmps;
                   
                            SSactc = [SSactc; 0 0 0];
                
                            ompt = ompt+jmpt.*SSactc(tmpt(:),:);
                            ompc = ompc+jmpc.*SSactc(tmpc(:),:);
                            omps = omps+jmps.*SSactc(tmps(:),:);           
                        end
                    end
                end
                
                curblob=st{curfig}.curblob;
                blobSeq=[curblob:numel(st{curfig}.vols{i}.blobs) 1:curblob-1];
                for j=blobSeq %1:numel(st{curfig}.vols{i}.blobs)
                    actc = st{curfig}.vols{i}.blobs{j}.colour.cmap;
                    actp = st{curfig}.vols{i}.blobs{j}.colour.prop;
                    % get max for blob image
                    if isfield(st{curfig}.vols{i}.blobs{j},'max')
                        cmx = st{curfig}.vols{i}.blobs{j}.max;
                    else
                        cmx = max([eps maxval(st{curfig}.vols{i}.blobs{j}.vol)]);
                    end
                    if isfield(st{curfig}.vols{i}.blobs{j},'min')
                        cmn = st{curfig}.vols{i}.blobs{j}.min;
                    else
                        cmn = -cmx;
                    end
                    
                    % get blob data
                    vol  = st{curfig}.vols{i}.blobs{j}.vol;
                    M    = st{curfig}.Space\st{curfig}.vols{i}.premul*st{curfig}.vols{i}.blobs{j}.mat;
                
                    %Revised by YAN Chao-Gan, 130720. Could also work with Data has been read into memory other than only depending on the file.
                    if ~isfield(vol,'Data')
                        tmpt = spm_slice_vol(vol,inv(TM0*M),TD,[0 NaN])';
                        tmpc = spm_slice_vol(vol,inv(CM0*M),CD,[0 NaN])';
                        tmps = spm_slice_vol(vol,inv(SM0*M),SD,[0 NaN])';
                    else   %Revised by YAN Chao-Gan, 130720. Could also work with Data has been read into memory other than only depending on the file.
                        tmpt = spm_slice_vol(vol.Data,inv(TM0*M),TD,[0 NaN])';
                        tmpc = spm_slice_vol(vol.Data,inv(CM0*M),CD,[0 NaN])';
                        tmps = spm_slice_vol(vol.Data,inv(SM0*M),SD,[0 NaN])';
                    end
                
                
                    % actimg scaled round 0, black NaNs
                    topc = size(actc,1)+1;
                    tmpt = scaletocmap(tmpt,cmn,cmx,actc,topc);
                    tmpc = scaletocmap(tmpc,cmn,cmx,actc,topc);
                    tmps = scaletocmap(tmps,cmn,cmx,actc,topc);
              
                    %Overlay Transparent Weight
                    jmpt = actp*repmat((tmpt(:)~=topc),1,3);
                    jmpc = actp*repmat((tmpc(:)~=topc),1,3);
                    jmps = actp*repmat((tmps(:)~=topc),1,3);
                    
                    %Except Underlay
                    umpt = umpt-jmpt;
                    umpc = umpc-jmpc;
                    umps = umps-jmps;
                    
                    %Negtive Weight Recoup
                    nmpt = umpt.*(umpt<0);
                    nmpc = umpc.*(umpc<0);
                    nmps = umps.*(umps<0);
                    
                    %Modified Overlay Transparent Weight
                    jmpt=jmpt+nmpt;
                    jmpc=jmpc+nmpc;
                    jmps=jmps+nmps;
                    
                    actc = [actc; 0 0 0];
                
                    ompt = ompt+jmpt.*actc(tmpt(:),:);
                    ompc = ompc+jmpc.*actc(tmpc(:),:);
                    omps = omps+jmps.*actc(tmps(:),:);
                    
                    umpt(umpt<0) = 0;
                    umpc(umpc<0) = 0;
                    umps(umps<0) = 0;
                end      
%                 if isfield(handles, 'DPABI_fig') && IsPreview
%                     SSactc=AStruct.Template.CMap;
%                     M = st{curfig}.Space\st{curfig}.vols{i}.premul*AStruct.Template.mat;
%                     
%                     
%                     tmpt = spm_slice_vol(Mask,inv(TM0*M),TD,[0 NaN])';
%                     tmpc = spm_slice_vol(Mask,inv(CM0*M),CD,[0 NaN])';
%                     tmps = spm_slice_vol(Mask,inv(SM0*M),SD,[0 NaN])';
%                     
%                     SStopc = size(SSactc,1)+1;
%                     tmpt = scaletocmap(tmpt, 0, 1, SSactc, SStopc);
%                     tmpc = scaletocmap(tmpc, 0, 1, SSactc, SStopc);
%                     tmps = scaletocmap(tmps, 0, 1, SSactc, SStopc);
%                     
%                     umpt = umpt+repmat((tmpt(:)~=SStopc),1,3);
%                     umpc = umpc+repmat((tmpc(:)~=SStopc),1,3);
%                     umps = umps+repmat((tmps(:)~=SStopc),1,3);
%                 
%                     SSactc = [SSactc; 0 0 0];
%                 
%                     ompt = ompt+SSactc(tmpt(:),:).*mmgt;
%                     ompc = ompc+SSactc(tmpc(:),:).*mmgc;
%                     omps = omps+SSactc(tmps(:),:).*mmgs;
%                 end
%                 
%                 umpt=umpt==0;
%                 umpc=umpc==0;
%                 umps=umps==0;
                
                
                % combine gray and blob data to
                % truecolour
                
                %Revised by Sandy, 20140104. 
                imgt = reshape(ompt+gryc(imgt(:),:).*umpt, [size(imgt) 3]);
                imgc = reshape(ompc+gryc(imgc(:),:).*umpc, [size(imgc) 3]);
                imgs = reshape(omps+gryc(imgs(:),:).*umps, [size(imgs) 3]);
                
                %Revised by YAN Chao-Gan, 130609. The unwanted voxels will keep full color of underlay.
%                 imgt = reshape(ompt*actp+ ...
%                     gryc(imgt(:),:)*(1-actp) + umpt.*gryc(imgt(:),:)*(actp)  , ...
%                     [size(imgt) 3]);
%                 imgc = reshape(ompc*actp+ ...
%                     gryc(imgc(:),:)*(1-actp) + umpc.*gryc(imgc(:),:)*(actp)  , ...
%                     [size(imgc) 3]);
%                 imgs = reshape(omps*actp+ ...
%                     gryc(imgs(:),:)*(1-actp) + umps.*gryc(imgs(:),:)*(actp)  , ...
%                     [size(imgs) 3]);
%                 imgt = reshape(actc(tmpt(:),:)*actp+ ...
%                     gryc(imgt(:),:)*(1-actp), ...
%                     [size(imgt) 3]);
%                 imgc = reshape(actc(tmpc(:),:)*actp+ ...
%                     gryc(imgc(:),:)*(1-actp), ...
%                     [size(imgc) 3]);
%                 imgs = reshape(actc(tmps(:),:)*actp+ ...
%                     gryc(imgs(:),:)*(1-actp), ...
%                     [size(imgs) 3]);
                %Revising finished, YAN Chao-Gan, 130609.
            else
                % Add full colour blobs - several sets at once
                scal  = 1/(mx-mn);
                dcoff = -mn*scal;
                
                wt = zeros(size(imgt));
                wc = zeros(size(imgc));
                ws = zeros(size(imgs));
                
                imgt  = repmat(imgt*scal+dcoff,[1,1,3]);
                imgc  = repmat(imgc*scal+dcoff,[1,1,3]);
                imgs  = repmat(imgs*scal+dcoff,[1,1,3]);
                
                cimgt = zeros(size(imgt));
                cimgc = zeros(size(imgc));
                cimgs = zeros(size(imgs));
                
                colour = zeros(numel(st{curfig}.vols{i}.blobs),3);
                for k=1:numel(st{curfig}.vols{i}.blobs) % get colours of all images first
                    if isfield(st{curfig}.vols{i}.blobs{k},'colour')
                        colour(k,:) = reshape(st{curfig}.vols{i}.blobs{k}.colour, [1 3]);
                    else
                        colour(k,:) = [1 0 0];
                    end
                end
                %colour = colour/max(sum(colour));
                
                for k=1:numel(st{curfig}.vols{i}.blobs)
                    if isfield(st{curfig}.vols{i}.blobs{k},'max')
                        mx = st{curfig}.vols{i}.blobs{k}.max;
                    else
                        mx = max([eps max(st{curfig}.vols{i}.blobs{k}.vol(:))]);
                        st{curfig}.vols{i}.blobs{k}.max = mx;
                    end
                    if isfield(st{curfig}.vols{i}.blobs{k},'min')
                        mn = st{curfig}.vols{i}.blobs{k}.min;
                    else
                        mn = min([0 min(st{curfig}.vols{i}.blobs{k}.vol(:))]);
                        st{curfig}.vols{i}.blobs{k}.min = mn;
                    end
                    
                    vol  = st{curfig}.vols{i}.blobs{k}.vol;
                    M    = st{curfig}.Space\st{curfig}.vols{i}.premul*st{curfig}.vols{i}.blobs{k}.mat;
                    
                    
                    if ~isfield(vol,'Data')
                        tmpt = spm_slice_vol(vol,inv(TM0*M),TD,[0 NaN])';
                        tmpc = spm_slice_vol(vol,inv(CM0*M),CD,[0 NaN])';
                        tmps = spm_slice_vol(vol,inv(SM0*M),SD,[0 NaN])';
                    else   %Revised by YAN Chao-Gan, 130720. Could also work with Data has been read into memory other than only depending on the file.
                        tmpt = spm_slice_vol(vol.Data,inv(TM0*M),TD,[0 NaN])';
                        tmpc = spm_slice_vol(vol.Data,inv(CM0*M),CD,[0 NaN])';
                        tmps = spm_slice_vol(vol.Data,inv(SM0*M),SD,[0 NaN])';
                    end
                    

                    % check min/max of sampled image
                    % against mn/mx as given in st
                    tmpt(tmpt(:)<mn) = mn;
                    tmpc(tmpc(:)<mn) = mn;
                    tmps(tmps(:)<mn) = mn;
                    tmpt(tmpt(:)>mx) = mx;
                    tmpc(tmpc(:)>mx) = mx;
                    tmps(tmps(:)>mx) = mx;
                    tmpt = (tmpt-mn)/(mx-mn);
                    tmpc = (tmpc-mn)/(mx-mn);
                    tmps = (tmps-mn)/(mx-mn);
                    tmpt(~isfinite(tmpt)) = 0;
                    tmpc(~isfinite(tmpc)) = 0;
                    tmps(~isfinite(tmps)) = 0;
                    
                    cimgt = cimgt + cat(3,tmpt*colour(k,1),tmpt*colour(k,2),tmpt*colour(k,3));
                    cimgc = cimgc + cat(3,tmpc*colour(k,1),tmpc*colour(k,2),tmpc*colour(k,3));
                    cimgs = cimgs + cat(3,tmps*colour(k,1),tmps*colour(k,2),tmps*colour(k,3));
                    
                    wt = wt + tmpt;
                    wc = wc + tmpc;
                    ws = ws + tmps;
                    cdata=permute(shiftdim((1/64:1/64:1)'* ...
                        colour(k,:),-1),[2 1 3]);
                    redraw_colourbar(i,k,[mn mx],cdata);
                end
                
                imgt = repmat(1-wt,[1 1 3]).*imgt+cimgt;
                imgc = repmat(1-wc,[1 1 3]).*imgc+cimgc;
                imgs = repmat(1-ws,[1 1 3]).*imgs+cimgs;
                
                imgt(imgt<0)=0; imgt(imgt>1)=1;
                imgc(imgc<0)=0; imgc(imgc>1)=1;
                imgs(imgs<0)=0; imgs(imgs>1)=1;
            end
        else
            if isfield(handles, 'DPABI_fig')
                IsPreview=get(handles.PreviewBtn, 'Value');
                if IsPreview
                    AtlasIdx=get(handles.StructuralPopup, 'Value');
                    AStruct=st{curfig}.AtlasInfo{AtlasIdx};
                    vol=AStruct.Template;
                    M=st{curfig}.Space\st{curfig}.vols{i}.premul*vol.mat;

                    AMat=AStruct.Template.mat;
                    APos=round(inv(AMat)*[X;Y;Z;1]);
                    AI=APos(1);
                    AJ=APos(2);
                    AK=APos(3);
            
                    AIndex=AStruct.Template.Data(AI, AJ, AK);
                                
                    if AIndex
                        Mask=double(AStruct.Template.Data==AIndex);
                        
                        tmpt = spm_slice_vol(Mask,inv(TM0*M),TD,[0 NaN])';
                        tmpc = spm_slice_vol(Mask,inv(CM0*M),CD,[0 NaN])';
                        tmps = spm_slice_vol(Mask,inv(SM0*M),SD,[0 NaN])';
                    else
                        IsPreview=0;
                    end
                end
            else
                IsPreview=0;
            end
            
            if IsPreview
                % colourmaps
                gryc = (0:63)'*ones(1,3)/63;
                
                % scale grayscale image, not isfinite -> black
                imgt = scaletocmap(imgt,mn,mx,gryc,65);
                imgc = scaletocmap(imgc,mn,mx,gryc,65);
                imgs = scaletocmap(imgs,mn,mx,gryc,65);
                
                gryc = [gryc; 0 0 0];
                
                SSactc=AStruct.Template.CMap;
                SSactp=0.8;
                    
                SStopc = size(SSactc,1)+1;
                
                tmpt = scaletocmap(tmpt, 0, 1, SSactc, SStopc);
                tmpc = scaletocmap(tmpc, 0, 1, SSactc, SStopc);
                tmps = scaletocmap(tmps, 0, 1, SSactc, SStopc);
                   
                umpt = repmat((tmpt(:)~=SStopc),1,3);
                umpc = repmat((tmpc(:)~=SStopc),1,3);
                umps = repmat((tmps(:)~=SStopc),1,3);
                
                SSactc = [SSactc; 0 0 0];
                
                ompt = SSactc(tmpt(:),:);
                ompc = SSactc(tmpc(:),:);
                omps = SSactc(tmps(:),:);
                
                umpt=umpt==0;
                umpc=umpc==0;
                umps=umps==0;
                
                imgt = reshape(ompt*SSactp+ ...
                    gryc(imgt(:),:)*(1-SSactp) + umpt.*gryc(imgt(:),:)*(SSactp)  , ...
                    [size(imgt) 3]);
                imgc = reshape(ompc*SSactp+ ...
                    gryc(imgc(:),:)*(1-SSactp) + umpc.*gryc(imgc(:),:)*(SSactp)  , ...
                    [size(imgc) 3]);
                imgs = reshape(omps*SSactp+ ...
                    gryc(imgs(:),:)*(1-SSactp) + umps.*gryc(imgs(:),:)*(SSactp)  , ...
                    [size(imgs) 3]);
            else
                scal = 64/(mx-mn);
                dcoff = -mn*scal;
                imgt = imgt*scal+dcoff;
                imgc = imgc*scal+dcoff;
                imgs = imgs*scal+dcoff;
            end
        end
        
        set(st{curfig}.vols{i}.ax{1}.d,'HitTest','off', 'Cdata',imgt);
        set(st{curfig}.vols{i}.ax{1}.lx,'HitTest','off',...
            'Xdata',[0 TD(1)]+0.5,'Ydata',[1 1]*(cent(2)-bb(1,2)+1));
        set(st{curfig}.vols{i}.ax{1}.ly,'HitTest','off',...
            'Ydata',[0 TD(2)]+0.5,'Xdata',[1 1]*(cent(1)-bb(1,1)+1));
        
        set(st{curfig}.vols{i}.ax{2}.d,'HitTest','off', 'Cdata',imgc);
        set(st{curfig}.vols{i}.ax{2}.lx,'HitTest','off',...
            'Xdata',[0 CD(1)]+0.5,'Ydata',[1 1]*(cent(3)-bb(1,3)+1));
        set(st{curfig}.vols{i}.ax{2}.ly,'HitTest','off',...
            'Ydata',[0 CD(2)]+0.5,'Xdata',[1 1]*(cent(1)-bb(1,1)+1));
        
        set(st{curfig}.vols{i}.ax{3}.d,'HitTest','off','Cdata',imgs);
        if st{curfig}.mode ==0
            set(st{curfig}.vols{i}.ax{3}.lx,'HitTest','off',...
                'Xdata',[0 SD(1)]+0.5,'Ydata',[1 1]*(cent(2)-bb(1,2)+1));
            set(st{curfig}.vols{i}.ax{3}.ly,'HitTest','off',...
                'Ydata',[0 SD(2)]+0.5,'Xdata',[1 1]*(cent(3)-bb(1,3)+1));
        else
            set(st{curfig}.vols{i}.ax{3}.lx,'HitTest','off',...
                'Xdata',[0 SD(1)]+0.5,'Ydata',[1 1]*(cent(3)-bb(1,3)+1));
            set(st{curfig}.vols{i}.ax{3}.ly,'HitTest','off',...
                'Ydata',[0 SD(2)]+0.5,'Xdata',[1 1]*(bb(2,2)+1-cent(2)));
        end
        
        if ~isempty(st{curfig}.plugins) % process any addons
            for k = 1:numel(st{curfig}.plugins)
                if isfield(st{curfig}.vols{i},st{curfig}.plugins{k})
                    feval(['w_spm_ov_', st{curfig}.plugins{k}], ...
                        'redraw', i, TM0, TD, CM0, CD, SM0, SD);
                end
            end
        end
    end
end
drawnow;


%==========================================================================
% function redraw_all
%==========================================================================
function redraw_all(varargin)
if ~isempty(varargin)
    redraw(varargin{1})
else
    curfig=GetCurFig;
    redraw(curfig);
end


%==========================================================================
% function redraw_colourbar(vh,bh,interval,cdata)
%==========================================================================
function redraw_colourbar(vh,bh,interval,cdata,curfig)
global st
if nargin<5
    curfig=GetCurFig;
end
if isfield(st{curfig}.vols{vh}.blobs{bh},'cbar')
    if st{curfig}.mode == 0
        axpos = get(st{curfig}.vols{vh}.ax{2}.ax,'Position');
    else
        axpos = get(st{curfig}.vols{vh}.ax{1}.ax,'Position');
    end
    % only scale cdata if we have out-of-range truecolour values
    if ndims(cdata)==3 && max(cdata(:))>1
        cdata=cdata./max(cdata(:));
    end
    image([0 1],interval,cdata,'Parent',st{curfig}.vols{vh}.blobs{bh}.cbar);
    handles=guidata(st{curfig}.fig);
    if ~isfield(handles, 'DPABI_fig')
        set(st{curfig}.vols{vh}.blobs{bh}.cbar, ...
            'Position',[(axpos(1)+axpos(3)+0.05+(bh-1)*.1)...
            (axpos(2)+0.005) 0.05 (axpos(4)-0.01)],...
            'YDir','normal','XTickLabel',[],'XTick',[]);
    else
        set(st{curfig}.vols{vh}.blobs{bh}.cbar, ...
            'YDir','normal','XTickLabel',[],'XTick',[]);
    end
    if isfield(st{curfig}.vols{vh}.blobs{bh},'name')
        ylabel(st{curfig}.vols{vh}.blobs{bh}.name,'parent',st{curfig}.vols{vh}.blobs{bh}.cbar);
    end
end


%==========================================================================
% function centre = findcent
%==========================================================================
function centre = findcent
global st
curfig=GetCurFig;
obj    = get(st{curfig}.fig,'CurrentObject');
centre = [];
cent   = [];
cp     = [];
for i=valid_handles
    for j=1:3
        if ~isempty(obj)
            if (st{curfig}.vols{i}.ax{j}.ax == obj),
                cp = get(obj,'CurrentPoint');
            end
        end
        if ~isempty(cp)
            cp   = cp(1,1:2);
            is   = inv(st{curfig}.Space);
            cent = is(1:3,1:3)*st{curfig}.centre(:) + is(1:3,4);
            switch j
                case 1
                    cent([1 2])=[cp(1)+st{curfig}.bb(1,1)-1 cp(2)+st{curfig}.bb(1,2)-1];
                case 2
                    cent([1 3])=[cp(1)+st{curfig}.bb(1,1)-1 cp(2)+st{curfig}.bb(1,3)-1];
                case 3
                    if st{curfig}.mode ==0
                        cent([3 2])=[cp(1)+st{curfig}.bb(1,3)-1 cp(2)+st{curfig}.bb(1,2)-1];
                    else
                        cent([2 3])=[st{curfig}.bb(2,2)+1-cp(1) cp(2)+st{curfig}.bb(1,3)-1];
                    end
            end
            break;
        end
    end
    if ~isempty(cent), break; end
end
if ~isempty(cent), centre = st{curfig}.Space(1:3,1:3)*cent(:) + st{curfig}.Space(1:3,4); end


%==========================================================================
% function handles = valid_handles(handles)
%==========================================================================
function handles = valid_handles(handles)
global st
curfig=GetCurFig;
if ~nargin, handles = 1:max_img; end
if isempty(st{curfig}) || ~isfield(st{curfig},'vols')
    handles = [];
else
    handles = handles(:)';
    handles = handles(handles<=max_img & handles>=1 & ~rem(handles,1));
    for h=handles
        if isempty(st{curfig}.vols{h}), handles(handles==h)=[]; end
    end
end


%==========================================================================
% function reset_st
%==========================================================================
function reset_st
global st
curfig=GetCurFig;
fig = spm_figure('FindWin','Graphics');
bb  = []; %[ [-78 78]' [-112 76]' [-50 85]' ];
st  = struct('n', 0, 'vols',{cell(max_img,1)}, 'bb',bb, 'Space',eye(4), ...
             'centre',[0 0 0], 'callback',';', 'xhairs',1, 'hld',1, ...
             'fig',fig, 'mode',1, 'plugins',{{}}, 'snap',[]);

xTB = spm('TBs');
if ~isempty(xTB)
    pluginbase = {spm('Dir') xTB.dir};
else
    pluginbase = {spm('Dir')};
end
for k = 1:numel(pluginbase)
    pluginpath = fullfile(pluginbase{k},'w_spm_orthviews');
    if isdir(pluginpath)
        pluginfiles = dir(fullfile(pluginpath,'spm_ov_*.m'));
        if ~isempty(pluginfiles)
            if ~isdeployed, addpath(pluginpath); end
            for l = 1:numel(pluginfiles)
                pluginname = spm_file(pluginfiles(l).name,'basename');
                st{curfig}.plugins{end+1} = strrep(pluginname, 'spm_ov_','');
            end
        end
    end
end


%==========================================================================
% function img = scaletocmap(inpimg,mn,mx,cmap,miscol)
%==========================================================================
function img = scaletocmap(inpimg,mn,mx,cmap,miscol)
if nargin < 5, miscol=1; end
cml = size(cmap,1);
scf = (cml-1)/(mx-mn);
img = round((inpimg-mn)*scf)+1;
img(img<1)   = 1; 
img(img>cml) = cml;
img(inpimg==0) = miscol; %Added by YAN Chao-Gan 130609, mask out the 0 voxels.
img(~isfinite(img)) = miscol;


%==========================================================================
% function cmap = getcmap(acmapname)
%==========================================================================
function cmap = getcmap(acmapname)
% get colormap of name acmapname
if ~isempty(acmapname)
    cmap = evalin('base',acmapname,'[]');
    if isempty(cmap) % not a matrix, is .mat file?
        acmat = spm_file(acmapname, 'ext','.mat');
        if exist(acmat, 'file')
            s    = struct2cell(load(acmat));
            cmap = s{1};
        end
    end
end
if size(cmap, 2)~=3
    warning('Colormap was not an N by 3 matrix')
    cmap = [];
end


%==========================================================================
% function item_parent = addcontext(volhandle)
%==========================================================================
function item_parent = addcontext(volhandle)
global st
curfig=GetCurFig;
% create context menu
set(0,'CurrentFigure',st{curfig}.fig);
% contextmenu
item_parent = uicontextmenu;

% contextsubmenu 0
item00 = uimenu(item_parent, 'Label','unknown image', 'UserData','filename');
w_spm_orthviews('context_menu','image_info',item00,volhandle);
item0a = uimenu(item_parent, 'UserData','pos_mm', 'Separator','on', ...
    'Callback','w_spm_orthviews(''context_menu'',''repos_mm'');');
item0b = uimenu(item_parent, 'UserData','pos_vx', ...
    'Callback','w_spm_orthviews(''context_menu'',''repos_vx'');');
item0c = uimenu(item_parent, 'UserData','v_value');

% contextsubmenu 1
item1    = uimenu(item_parent,'Label','Zoom', 'Separator','on');
[zl, rl] = w_spm_orthviews('ZoomMenu');
for cz = numel(zl):-1:1
    if isinf(zl(cz))
        czlabel = 'Full Volume';
    elseif isnan(zl(cz))
        czlabel = 'BBox, this image > ...';
    elseif zl(cz) == 0
        czlabel = 'BBox, this image nonzero';
    else
        czlabel = sprintf('%dx%d mm', 2*zl(cz), 2*zl(cz));
    end
    item1_x = uimenu(item1, 'Label',czlabel,...
        'Callback', sprintf(...
        'w_spm_orthviews(''context_menu'',''zoom'',%d,%d)',zl(cz),rl(cz)));
    if isinf(zl(cz)) % default display is Full Volume
        set(item1_x, 'Checked','on');
    end
end

% contextsubmenu 2
checked   = {'off','off'};
checked{st{curfig}.xhairs+1} = 'on';
item2     = uimenu(item_parent,'Label','Crosshairs');
item2_1   = uimenu(item2,      'Label','on',  'Callback','w_spm_orthviews(''context_menu'',''Xhair'',''on'');','Checked',checked{2});
item2_2   = uimenu(item2,      'Label','off', 'Callback','w_spm_orthviews(''context_menu'',''Xhair'',''off'');','Checked',checked{1});

% contextsubmenu 3
if st{curfig}.Space == eye(4)
    checked = {'off', 'on'};
else
    checked = {'on', 'off'};
end
item3     = uimenu(item_parent,'Label','Orientation');
item3_1   = uimenu(item3,      'Label','World space', 'Callback','w_spm_orthviews(''context_menu'',''orientation'',3);','Checked',checked{2});
item3_2   = uimenu(item3,      'Label','Voxel space (1st image)', 'Callback','w_spm_orthviews(''context_menu'',''orientation'',2);','Checked',checked{1});
item3_3   = uimenu(item3,      'Label','Voxel space (this image)', 'Callback','w_spm_orthviews(''context_menu'',''orientation'',1);','Checked','off');

% contextsubmenu 3
if isempty(st{curfig}.snap)
    checked = {'off', 'on'};
else
    checked = {'on', 'off'};
end
item3     = uimenu(item_parent,'Label','Snap to Grid');
item3_1   = uimenu(item3,      'Label','Don''t snap', 'Callback','w_spm_orthviews(''context_menu'',''snap'',3);','Checked',checked{2});
item3_2   = uimenu(item3,      'Label','Snap to 1st image', 'Callback','w_spm_orthviews(''context_menu'',''snap'',2);','Checked',checked{1});
item3_3   = uimenu(item3,      'Label','Snap to this image', 'Callback','w_spm_orthviews(''context_menu'',''snap'',1);','Checked','off');

% contextsubmenu 4
if st{curfig}.hld == 0
    checked = {'off', 'off', 'on'};
elseif st{curfig}.hld > 0
    checked = {'off', 'on', 'off'};
else
    checked = {'on', 'off', 'off'};
end
item4     = uimenu(item_parent,'Label','Interpolation');
item4_1   = uimenu(item4,      'Label','NN',    'Callback','w_spm_orthviews(''context_menu'',''interpolation'',3);', 'Checked',checked{3});
item4_2   = uimenu(item4,      'Label','Trilin', 'Callback','w_spm_orthviews(''context_menu'',''interpolation'',2);','Checked',checked{2});
item4_3   = uimenu(item4,      'Label','Sinc',  'Callback','w_spm_orthviews(''context_menu'',''interpolation'',1);','Checked',checked{1});

% contextsubmenu 5
% item5     = uimenu(item_parent,'Label','Position', 'Callback','w_spm_orthviews(''context_menu'',''position'');');

% contextsubmenu 6
item6       = uimenu(item_parent,'Label','Image','Separator','on');
item6_1     = uimenu(item6,      'Label','Window');
item6_1_1   = uimenu(item6_1,    'Label','local');
item6_1_1_1 = uimenu(item6_1_1,  'Label','auto', 'Callback','w_spm_orthviews(''context_menu'',''window'',2);');
item6_1_1_2 = uimenu(item6_1_1,  'Label','manual', 'Callback','w_spm_orthviews(''context_menu'',''window'',1);');
item6_1_1_3 = uimenu(item6_1_1,  'Label','percentiles', 'Callback','w_spm_orthviews(''context_menu'',''window'',3);');
item6_1_2   = uimenu(item6_1,    'Label','global');
item6_1_2_1 = uimenu(item6_1_2,  'Label','auto', 'Callback','w_spm_orthviews(''context_menu'',''window_gl'',2);');
item6_1_2_2 = uimenu(item6_1_2,  'Label','manual', 'Callback','w_spm_orthviews(''context_menu'',''window_gl'',1);');
if license('test','image_toolbox') == 1
    offon = {'off', 'on'};
    checked = offon(strcmp(st{curfig}.vols{volhandle}.mapping, ...
        {'linear', 'histeq', 'loghisteq', 'quadhisteq'})+1);
    item6_2     = uimenu(item6,      'Label','Intensity mapping');
    item6_2_1   = uimenu(item6_2,    'Label','local');
    item6_2_1_1 = uimenu(item6_2_1,  'Label','Linear', 'Checked',checked{1}, ...
        'Callback','w_spm_orthviews(''context_menu'',''mapping'',''linear'');');
    item6_2_1_2 = uimenu(item6_2_1,  'Label','Equalised histogram', 'Checked',checked{2}, ...
        'Callback','w_spm_orthviews(''context_menu'',''mapping'',''histeq'');');
    item6_2_1_3 = uimenu(item6_2_1,  'Label','Equalised log-histogram', 'Checked',checked{3}, ...
        'Callback','w_spm_orthviews(''context_menu'',''mapping'',''loghisteq'');');
    item6_2_1_4 = uimenu(item6_2_1,  'Label','Equalised squared-histogram', 'Checked',checked{4}, ...
        'Callback','w_spm_orthviews(''context_menu'',''mapping'',''quadhisteq'');');
    item6_2_2   = uimenu(item6_2,    'Label','global');
    item6_2_2_1 = uimenu(item6_2_2,  'Label','Linear', 'Checked',checked{1}, ...
        'Callback','w_spm_orthviews(''context_menu'',''mapping_gl'',''linear'');');
    item6_2_2_2 = uimenu(item6_2_2,  'Label','Equalised histogram', 'Checked',checked{2}, ...
        'Callback','w_spm_orthviews(''context_menu'',''mapping_gl'',''histeq'');');
    item6_2_2_3 = uimenu(item6_2_2,  'Label','Equalised log-histogram', 'Checked',checked{3}, ...
        'Callback','w_spm_orthviews(''context_menu'',''mapping_gl'',''loghisteq'');');
    item6_2_2_4 = uimenu(item6_2_2,  'Label','Equalised squared-histogram', 'Checked',checked{4}, ...
        'Callback','w_spm_orthviews(''context_menu'',''mapping_gl'',''quadhisteq'');');
end

% contextsubmenu 7
item7     = uimenu(item_parent,'Label','Blobs');
item7_1   = uimenu(item7,      'Label','Add blobs');
item7_1_1 = uimenu(item7_1,    'Label','local',  'Callback','w_spm_orthviews(''context_menu'',''add_blobs'',2);');
item7_1_2 = uimenu(item7_1,    'Label','global', 'Callback','w_spm_orthviews(''context_menu'',''add_blobs'',1);');
item7_2   = uimenu(item7,      'Label','Add image');
item7_2_1 = uimenu(item7_2,    'Label','local',  'Callback','w_spm_orthviews(''context_menu'',''add_image'',2);');
item7_2_2 = uimenu(item7_2,    'Label','global', 'Callback','w_spm_orthviews(''context_menu'',''add_image'',1);');
item7_3   = uimenu(item7,      'Label','Add colored blobs','Separator','on');
item7_3_1 = uimenu(item7_3,    'Label','local',  'Callback','w_spm_orthviews(''context_menu'',''add_c_blobs'',2);');
item7_3_2 = uimenu(item7_3,    'Label','global', 'Callback','w_spm_orthviews(''context_menu'',''add_c_blobs'',1);');
item7_4   = uimenu(item7,      'Label','Add colored image');
item7_4_1 = uimenu(item7_4,    'Label','local',  'Callback','w_spm_orthviews(''context_menu'',''add_c_image'',2);');
item7_4_2 = uimenu(item7_4,    'Label','global', 'Callback','w_spm_orthviews(''context_menu'',''add_c_image'',1);');
item7_5   = uimenu(item7,      'Label','Remove blobs',        'Visible','off','Separator','on');
item7_6   = uimenu(item7,      'Label','Remove colored blobs','Visible','off');
item7_6_1 = uimenu(item7_6,    'Label','local', 'Visible','on');
item7_6_2 = uimenu(item7_6,    'Label','global','Visible','on');
item7_7   = uimenu(item7,      'Label','Set blobs max', 'Visible','off');

for i=1:3
    set(st{curfig}.vols{volhandle}.ax{i}.ax,'UIcontextmenu',item_parent);
    st{curfig}.vols{volhandle}.ax{i}.cm = item_parent;
end

% process any plugins
for k = 1:numel(st{curfig}.plugins)
    feval(['w_spm_ov_', st{curfig}.plugins{k}],'context_menu',volhandle,item_parent);
    if k==1
        h = get(item_parent,'Children');
        set(h(1),'Separator','on'); 
    end
end


%==========================================================================
% function addcontexts(handles)
%==========================================================================
function addcontexts(handles)
for ii = valid_handles(handles)
    addcontext(ii);
end
w_spm_orthviews('reposition',w_spm_orthviews('pos'));


%==========================================================================
% function rmcontexts(handles)
%==========================================================================
function rmcontexts(handles)
global st
curfig=GetCurFig;
for ii = valid_handles(handles)
    for i=1:3
        set(st{curfig}.vols{ii}.ax{i}.ax,'UIcontextmenu',[]);
        try st{curfig}.vols{ii}.ax{i} = rmfield(st{curfig}.vols{ii}.ax{i},'cm'); end
    end
end


%==========================================================================
% function c_menu(varargin)
%==========================================================================
function c_menu(varargin)
global st
curfig=GetCurFig;

switch lower(varargin{1})
    case 'image_info'
        if nargin <3
            current_handle = get_current_handle;
        else
            current_handle = varargin{3};
        end
        if isfield(st{curfig}.vols{current_handle},'fname')
            [p,n,e,v] = spm_fileparts(st{curfig}.vols{current_handle}.fname);
            if isfield(st{curfig}.vols{current_handle},'n')
                v = sprintf(',%d',st{curfig}.vols{current_handle}.n);
            end
            set(varargin{2}, 'Label',[n e v]);
        end
        delete(get(varargin{2},'children'));
        if exist('p','var')
            item1 = uimenu(varargin{2}, 'Label', p);
        end
        if isfield(st{curfig}.vols{current_handle},'descrip')
            item2 = uimenu(varargin{2}, 'Label',...
                st{curfig}.vols{current_handle}.descrip);
        end
        dt = st{curfig}.vols{current_handle}.dt(1);
        item3 = uimenu(varargin{2}, 'Label', sprintf('Data type: %s', spm_type(dt)));
        str   = 'Intensity: varied';
        if size(st{curfig}.vols{current_handle}.pinfo,2) == 1
            if st{curfig}.vols{current_handle}.pinfo(2)
                str = sprintf('Intensity: Y = %g X + %g',...
                    st{curfig}.vols{current_handle}.pinfo(1:2)');
            else
                str = sprintf('Intensity: Y = %g X', st{curfig}.vols{current_handle}.pinfo(1)');
            end
        end
        item4  = uimenu(varargin{2}, 'Label',str);
        item5  = uimenu(varargin{2}, 'Label', 'Image dimensions', 'Separator','on');
        item51 = uimenu(varargin{2}, 'Label',...
            sprintf('%dx%dx%d', st{curfig}.vols{current_handle}.dim(1:3)));
        
        prms   = spm_imatrix(st{curfig}.vols{current_handle}.mat);
        item6  = uimenu(varargin{2}, 'Label', 'Voxel size', 'Separator','on');
        item61 = uimenu(varargin{2}, 'Label', sprintf('%.2f %.2f %.2f', prms(7:9)));
        
        O      = st{curfig}.vols{current_handle}.mat\[0 0 0 1]'; O=O(1:3)';
        item7  = uimenu(varargin{2}, 'Label', 'Origin', 'Separator','on');
        item71 = uimenu(varargin{2}, 'Label', sprintf('%.2f %.2f %.2f', O));
        
        R      = spm_matrix([0 0 0 prms(4:6)]);
        item8  = uimenu(varargin{2}, 'Label', 'Rotations', 'Separator','on');
        item81 = uimenu(varargin{2}, 'Label', sprintf('%.2f %.2f %.2f', R(1,1:3)));
        item82 = uimenu(varargin{2}, 'Label', sprintf('%.2f %.2f %.2f', R(2,1:3)));
        item83 = uimenu(varargin{2}, 'Label', sprintf('%.2f %.2f %.2f', R(3,1:3)));
        item9  = uimenu(varargin{2},...
            'Label','Specify other image...',...
            'Callback','w_spm_orthviews(''context_menu'',''swap_img'');',...
            'Separator','on');
        
    case 'repos_mm'
        oldpos_mm = w_spm_orthviews('pos');
        newpos_mm = spm_input('New Position (mm)','+1','r',sprintf('%.2f %.2f %.2f',oldpos_mm),3);
        w_spm_orthviews('reposition',newpos_mm);
        
    case 'repos_vx'
        current_handle = get_current_handle;
        oldpos_vx = w_spm_orthviews('pos', current_handle);
        newpos_vx = spm_input('New Position (voxels)','+1','r',sprintf('%.2f %.2f %.2f',oldpos_vx),3);
        newpos_mm = st{curfig}.vols{current_handle}.mat*[newpos_vx;1];
        w_spm_orthviews('reposition',newpos_mm(1:3));
        
    case 'zoom'
        zoom_all(varargin{2:end});
        bbox;
        redraw_all;
        
    case 'xhair'
        w_spm_orthviews('Xhairs',varargin{2});
        cm_handles = get_cm_handles;
        for i = 1:numel(cm_handles)
            z_handle = get(findobj(cm_handles(i),'label','Crosshairs'),'Children');
            set(z_handle,'Checked','off'); %reset check
            if strcmp(varargin{2},'off'), op = 1; else op = 2; end
            set(z_handle(op),'Checked','on');
        end
        
    case 'orientation'
        cm_handles = get_cm_handles;
        for i = 1:numel(cm_handles)
            z_handle = get(findobj(cm_handles(i),'label','Orientation'),'Children');
            set(z_handle,'Checked','off');
        end
        if varargin{2} == 3
            w_spm_orthviews('Space');
            for i = 1:numel(cm_handles),
                z_handle = findobj(cm_handles(i),'label','World space');
                set(z_handle,'Checked','on');
            end
        elseif varargin{2} == 2,
            w_spm_orthviews('Space',1);
            for i = 1:numel(cm_handles)
                z_handle = findobj(cm_handles(i),'label',...
                    'Voxel space (1st image)');
                set(z_handle,'Checked','on');
            end
        else
            w_spm_orthviews('Space',get_current_handle);
            z_handle = findobj(st{curfig}.vols{get_current_handle}.ax{1}.cm, ...
                'label','Voxel space (this image)');
            set(z_handle,'Checked','on');
            return;
        end
        
    case 'snap'
        cm_handles = get_cm_handles;
        for i = 1:numel(cm_handles)
            z_handle = get(findobj(cm_handles(i),'label','Snap to Grid'),'Children');
            set(z_handle,'Checked','off');
        end
        if varargin{2} == 3
            st{curfig}.snap = [];
        elseif varargin{2} == 2
            st{curfig}.snap = 1;
        else
            st{curfig}.snap = get_current_handle;
            z_handle = get(findobj(st{curfig}.vols{get_current_handle}.ax{1}.cm,'label','Snap to Grid'),'Children');
            set(z_handle(1),'Checked','on');
            return;
        end
        for i = 1:numel(cm_handles)
            z_handle = get(findobj(cm_handles(i),'label','Snap to Grid'),'Children');
            set(z_handle(varargin{2}),'Checked','on');
        end
        
    case 'interpolation'
        tmp        = [-4 1 0];
        st{curfig}.hld     = tmp(varargin{2});
        cm_handles = get_cm_handles;
        for i = 1:numel(cm_handles)
            z_handle = get(findobj(cm_handles(i),'label','Interpolation'),'Children');
            set(z_handle,'Checked','off');
            set(z_handle(varargin{2}),'Checked','on');
        end
        redraw_all;
        
    case 'window'
        current_handle = get_current_handle;
        if varargin{2} == 2
            w_spm_orthviews('window',current_handle);
        elseif varargin{2} == 3
            pc = spm_input('Percentiles', '+1', 'w', '3 97', 2, 100);
            wn = spm_summarise(st{curfig}.vols{current_handle}, 'all', ...
                @(X) spm_percentile(X, pc));
            w_spm_orthviews('window',current_handle,wn);
        else
            if isnumeric(st{curfig}.vols{current_handle}.window)
                defstr = sprintf('%.2f %.2f', st{curfig}.vols{current_handle}.window);
            else
                defstr = '';
            end
            [w,yp] = spm_input('Range','+1','e',defstr,[1 inf]);
            while numel(w) < 1 || numel(w) > 2
                uiwait(warndlg('Window must be one or two numbers','Wrong input size','modal'));
                [w,yp] = spm_input('Range',yp,'e',defstr,[1 inf]);
            end
            if numel(w) == 1
                w(2) = w(1)+eps;
            end
            w_spm_orthviews('window',current_handle,w);
        end
        
    case 'window_gl'
        if varargin{2} == 2
            for i = 1:numel(get_cm_handles)
                st{curfig}.vols{i}.window = 'auto';
            end
        else
            current_handle = get_current_handle;
            if isnumeric(st{curfig}.vols{current_handle}.window)
                defstr = sprintf('%d %d', st{curfig}.vols{current_handle}.window);
            else
                defstr = '';
            end
            [w,yp] = spm_input('Range','+1','e',defstr,[1 inf]);
            while numel(w) < 1 || numel(w) > 2
                uiwait(warndlg('Window must be one or two numbers','Wrong input size','modal'));
                [w,yp] = spm_input('Range',yp,'e',defstr,[1 inf]);
            end
            if numel(w) == 1
                w(2) = w(1)+eps;
            end
            for i = 1:numel(get_cm_handles)
                st{curfig}.vols{i}.window = w;
            end
        end
        redraw_all;
        
    case 'mapping'
        checked = strcmp(varargin{2}, ...
            {'linear', 'histeq', 'loghisteq', 'quadhisteq'});
        checked = checked(end:-1:1); % Handles are stored in inverse order
        current_handle = get_current_handle;
        cm_handles = get_cm_handles;
        st{curfig}.vols{current_handle}.mapping = varargin{2};
        z_handle = get(findobj(cm_handles(current_handle), ...
            'label','Intensity mapping'),'Children');
        for k = 1:numel(z_handle)
            c_handle = get(z_handle(k), 'Children');
            set(c_handle, 'checked', 'off');
            set(c_handle(checked), 'checked', 'on');
        end
        redraw_all;
        
    case 'mapping_gl'
        checked = strcmp(varargin{2}, ...
            {'linear', 'histeq', 'loghisteq', 'quadhisteq'});
        checked = checked(end:-1:1); % Handles are stored in inverse order
        cm_handles = get_cm_handles;
        for k = valid_handles
            st{curfig}.vols{k}.mapping = varargin{2};
            z_handle = get(findobj(cm_handles(k), ...
                'label','Intensity mapping'),'Children');
            for l = 1:numel(z_handle)
                c_handle = get(z_handle(l), 'Children');
                set(c_handle, 'checked', 'off');
                set(c_handle(checked), 'checked', 'on');
            end
        end
        redraw_all;
        
    case 'swap_img'
        current_handle = get_current_handle;
        newimg = spm_select(1,'image','select new image');
        if ~isempty(newimg)
            new_info = spm_vol(newimg);
            fn = fieldnames(new_info);
            for k=1:numel(fn)
                st{curfig}.vols{current_handle}.(fn{k}) = new_info.(fn{k});
            end
            w_spm_orthviews('context_menu','image_info',get(gcbo, 'parent'));
            redraw_all;
        end
        
    case 'add_blobs'
        % Add blobs to the image - in split colortable
        cm_handles = valid_handles;
        if varargin{2} == 2, cm_handles = get_current_handle; end
        spm_input('!DeleteInputObj');
        [SPM,xSPM] = spm_getSPM;
        if ~isempty(SPM)
            for i = 1:numel(cm_handles)
                addblobs(cm_handles(i),xSPM.XYZ,xSPM.Z,xSPM.M);
                % Add options for removing blobs
                c_handle = findobj(findobj(st{curfig}.vols{cm_handles(i)}.ax{1}.cm,'label','Blobs'),'Label','Remove blobs');
                set(c_handle,'Visible','on');
                delete(get(c_handle,'Children'));
                item7_3_1 = uimenu(c_handle,'Label','local','Callback','w_spm_orthviews(''context_menu'',''remove_blobs'',2);');
                if varargin{2} == 1,
                    item7_3_2 = uimenu(c_handle,'Label','global','Callback','w_spm_orthviews(''context_menu'',''remove_blobs'',1);');
                end
                % Add options for setting maxima for blobs
                c_handle = findobj(findobj(st{curfig}.vols{cm_handles(i)}.ax{1}.cm,'label','Blobs'),'Label','Set blobs max');
                set(c_handle,'Visible','on');
                delete(get(c_handle,'Children'));
                uimenu(c_handle,'Label','local','Callback','w_spm_orthviews(''context_menu'',''setblobsmax'',2);');
                if varargin{2} == 1
                    uimenu(c_handle,'Label','global','Callback','w_spm_orthviews(''context_menu'',''setblobsmax'',1);');
                end
            end
            redraw_all;
        end
        
    case 'remove_blobs'
        cm_handles = valid_handles;
        if varargin{2} == 2, cm_handles = get_current_handle; end
        for i = 1:numel(cm_handles)
            rmblobs(cm_handles(i));
            % Remove options for removing blobs
            c_handle = findobj(findobj(st{curfig}.vols{cm_handles(i)}.ax{1}.cm,'label','Blobs'),'Label','Remove blobs');
            delete(get(c_handle,'Children'));
            set(c_handle,'Visible','off');
            % Remove options for setting maxima for blobs
            c_handle = findobj(findobj(st{curfig}.vols{cm_handles(i)}.ax{1}.cm,'label','Blobs'),'Label','Set blobs max');
            set(c_handle,'Visible','off');
        end
        redraw_all;
        
    case 'add_image'
        % Add blobs to the image - in split colortable
        cm_handles = valid_handles;
        if varargin{2} == 2, cm_handles = get_current_handle; end
        spm_input('!DeleteInputObj');
        fname = spm_select(1,'image','select image');
        if ~isempty(fname)
            for i = 1:numel(cm_handles)
                addimage(cm_handles(i),fname);
                % Add options for removing blobs
                c_handle = findobj(findobj(st{curfig}.vols{cm_handles(i)}.ax{1}.cm,'label','Blobs'),'Label','Remove blobs');
                set(c_handle,'Visible','on');
                delete(get(c_handle,'Children'));
                item7_3_1 = uimenu(c_handle,'Label','local','Callback','w_spm_orthviews(''context_menu'',''remove_blobs'',2);');
                if varargin{2} == 1
                    item7_3_2 = uimenu(c_handle,'Label','global','Callback','w_spm_orthviews(''context_menu'',''remove_blobs'',1);');
                end
                % Add options for setting maxima for blobs
                c_handle = findobj(findobj(st{curfig}.vols{cm_handles(i)}.ax{1}.cm,'label','Blobs'),'Label','Set blobs max');
                set(c_handle,'Visible','on');
                delete(get(c_handle,'Children'));
                uimenu(c_handle,'Label','local','Callback','w_spm_orthviews(''context_menu'',''setblobsmax'',2);');
                if varargin{2} == 1
                    uimenu(c_handle,'Label','global','Callback','w_spm_orthviews(''context_menu'',''setblobsmax'',1);');
                end
            end
            redraw_all;
        end
        
    case 'add_c_blobs'
        % Add blobs to the image - in full colour
        cm_handles = valid_handles;
        if varargin{2} == 2, cm_handles = get_current_handle; end
        spm_input('!DeleteInputObj');
        [SPM,xSPM] = spm_getSPM;
        if ~isempty(SPM)
            c = spm_input('Colour','+1','m',...
                'Red blobs|Yellow blobs|Green blobs|Cyan blobs|Blue blobs|Magenta blobs',[1 2 3 4 5 6],1);
            colours = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1];
            c_names = {'red';'yellow';'green';'cyan';'blue';'magenta'};
            hlabel = sprintf('%s (%s)',xSPM.title,c_names{c});
            for i = 1:numel(cm_handles)
                addcolouredblobs(cm_handles(i),xSPM.XYZ,xSPM.Z,xSPM.M,colours(c,:),xSPM.title);
                addcolourbar(cm_handles(i),numel(st{curfig}.vols{cm_handles(i)}.blobs));
                c_handle    = findobj(findobj(st{curfig}.vols{cm_handles(i)}.ax{1}.cm,'label','Blobs'),'Label','Remove colored blobs');
                ch_c_handle = get(c_handle,'Children');
                set(c_handle,'Visible','on');
                %set(ch_c_handle,'Visible',on');
                item7_4_1   = uimenu(ch_c_handle(2),'Label',hlabel,'ForegroundColor',colours(c,:),...
                    'Callback','c = get(gcbo,''UserData'');w_spm_orthviews(''context_menu'',''remove_c_blobs'',2,c);',...
                    'UserData',c);
                if varargin{2} == 1
                    item7_4_2 = uimenu(ch_c_handle(1),'Label',hlabel,'ForegroundColor',colours(c,:),...
                        'Callback','c = get(gcbo,''UserData'');w_spm_orthviews(''context_menu'',''remove_c_blobs'',1,c);',...
                        'UserData',c);
                end
            end
            redraw_all;
        end
        
    case 'remove_c_blobs'
        cm_handles = valid_handles;
        if varargin{2} == 2, cm_handles = get_current_handle; end
        colours = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1];
        for i = 1:numel(cm_handles)
            if isfield(st{curfig}.vols{cm_handles(i)},'blobs')
                for j = 1:numel(st{curfig}.vols{cm_handles(i)}.blobs)
                    if all(st{curfig}.vols{cm_handles(i)}.blobs{j}.colour == colours(varargin{3},:));
                        if isfield(st{curfig}.vols{cm_handles(i)}.blobs{j},'cbar')
                            delete(st{curfig}.vols{cm_handles(i)}.blobs{j}.cbar);
                        end
                        st{curfig}.vols{cm_handles(i)}.blobs(j) = [];
                        break;
                    end
                end
                rm_c_menu = findobj(st{curfig}.vols{cm_handles(i)}.ax{1}.cm,'Label','Remove colored blobs');
                delete(gcbo);
                if isempty(st{curfig}.vols{cm_handles(i)}.blobs)
                    st{curfig}.vols{cm_handles(i)} = rmfield(st{curfig}.vols{cm_handles(i)},'blobs');
                    set(rm_c_menu, 'Visible', 'off');
                end
            end
        end
        redraw_all;
        
    case 'add_c_image'
        % Add truecolored image
        cm_handles = valid_handles;
        if varargin{2} == 2, cm_handles = get_current_handle; end
        spm_input('!DeleteInputObj');
        fname = spm_select([1 Inf],'image','select image(s)');
        for k = 1:size(fname,1)
            c = spm_input(sprintf('Image %d: Colour',k),'+1','m',...
                'Red blobs|Yellow blobs|Green blobs|Cyan blobs|Blue blobs|Magenta blobs',[1 2 3 4 5 6],1);
            colours = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1];
            c_names = {'red';'yellow';'green';'cyan';'blue';'magenta'};
            hlabel = sprintf('%s (%s)',fname(k,:),c_names{c});
            for i = 1:numel(cm_handles)
                addcolouredimage(cm_handles(i),fname(k,:),colours(c,:));
                addcolourbar(cm_handles(i),numel(st{curfig}.vols{cm_handles(i)}.blobs));
                c_handle    = findobj(findobj(st{curfig}.vols{cm_handles(i)}.ax{1}.cm,'label','Blobs'),'Label','Remove colored blobs');
                ch_c_handle = get(c_handle,'Children');
                set(c_handle,'Visible','on');
                %set(ch_c_handle,'Visible',on');
                item7_4_1 = uimenu(ch_c_handle(2),'Label',hlabel,'ForegroundColor',colours(c,:),...
                    'Callback','c = get(gcbo,''UserData'');w_spm_orthviews(''context_menu'',''remove_c_blobs'',2,c);','UserData',c);
                if varargin{2} == 1
                    item7_4_2 = uimenu(ch_c_handle(1),'Label',hlabel,'ForegroundColor',colours(c,:),...
                        'Callback','c = get(gcbo,''UserData'');w_spm_orthviews(''context_menu'',''remove_c_blobs'',1,c);',...
                        'UserData',c);
                end
            end
            redraw_all;
        end
        
    case 'setblobsmax'
        if varargin{2} == 1
            % global
            cm_handles = valid_handles;
            mx = -inf;
            for i = 1:numel(cm_handles)
                if ~isfield(st{curfig}.vols{cm_handles(i)}, 'blobs'), continue, end
                for j = 1:numel(st{curfig}.vols{cm_handles(i)}.blobs)
                    mx = max(mx, st{curfig}.vols{cm_handles(i)}.blobs{j}.max);
                end
            end
            mx = spm_input('Maximum value', '+1', 'r', mx, 1);
            for i = 1:numel(cm_handles)
                if ~isfield(st{curfig}.vols{cm_handles(i)}, 'blobs'), continue, end
                for j = 1:numel(st{curfig}.vols{cm_handles(i)}.blobs)
                    st{curfig}.vols{cm_handles(i)}.blobs{j}.max = mx;
                end
            end
        else
            % local (should handle coloured blobs, but not implemented yet)
            cm_handle = get_current_handle;
            colours = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1];
            if ~isfield(st{curfig}.vols{cm_handle}, 'blobs'), return, end
            for j = 1:numel(st{curfig}.vols{cm_handle}.blobs)
                if nargin < 4 || ...
                        all(st{curfig}.vols{cm_handle}.blobs{j}.colour == colours(varargin{3},:))
                    mx = st{curfig}.vols{cm_handle}.blobs{j}.max;
                    mx = spm_input('Maximum value', '+1', 'r', mx, 1);
                    st{curfig}.vols{cm_handle}.blobs{j}.max = mx;
                end
            end
        end
        redraw_all;
end


%==========================================================================
% function current_handle = get_current_handle
%==========================================================================
function current_handle = get_current_handle
cm_handle      = get(gca,'UIContextMenu');
cm_handles     = get_cm_handles;
current_handle = find(cm_handles==cm_handle);


%==========================================================================
% function cm_pos
%==========================================================================
function cm_pos(varargin)
global st
if ~isempty(varargin)
    curfig=varargin{1};
else
    curfig=GetCurFig;
end
for i = 1:numel(valid_handles)
    if isfield(st{curfig}.vols{i}.ax{1},'cm')
        set(findobj(st{curfig}.vols{i}.ax{1}.cm,'UserData','pos_mm'),...
            'Label',sprintf('mm:  %.1f %.1f %.1f',w_spm_orthviews('pos')));
        pos = w_spm_orthviews('pos',i);
        set(findobj(st{curfig}.vols{i}.ax{1}.cm,'UserData','pos_vx'),...
            'Label',sprintf('vx:  %.1f %.1f %.1f',pos));
        try
            %Fixed by Sandy 20140106
            if isfield(st{curfig}.vols{i}, 'Data')
                Y = spm_sample_vol(st{curfig}.vols{i}.Data,pos(1),pos(2),pos(3),st{curfig}.hld);
            else
                Y = spm_sample_vol(st{curfig}.vols{i},pos(1),pos(2),pos(3),st{curfig}.hld);
            end
        catch
            Y = NaN;
            fprintf('Cannot access file "%s".\n', st{curfig}.vols{i}.fname);
        end
        set(findobj(st{curfig}.vols{i}.ax{1}.cm,'UserData','v_value'),...
            'Label',sprintf('Y = %g',Y));
    end
end


%==========================================================================
% function cm_handles = get_cm_handles
%==========================================================================
function cm_handles = get_cm_handles
global st
curfig=GetCurFig;
cm_handles = [];
for i = valid_handles
    cm_handles = [cm_handles st{curfig}.vols{i}.ax{1}.cm];
end


%==========================================================================
% function zoom_all(zoom,res)
%==========================================================================
function zoom_all(zoom,res)
cm_handles = get_cm_handles;
zoom_op(zoom,res);
for i = 1:numel(cm_handles)
    z_handle = get(findobj(cm_handles(i),'label','Zoom'),'Children');
    set(z_handle,'Checked','off');
    if isinf(zoom)
        set(findobj(z_handle,'Label','Full Volume'),'Checked','on');
    elseif zoom > 0
        set(findobj(z_handle,'Label',sprintf('%dx%d mm', 2*zoom, 2*zoom)),'Checked','on');
    end % leave all unchecked if either bounding box option was chosen
end


%==========================================================================
% function m = max_img
%==========================================================================
function m = max_img
m = 1;

function curfig = GetCurFig(varargin)
if nargin<2
    curfig=gcf;
    curfig=w_Compatible2014bFig(curfig);
    if rem(curfig, 1)
        curfig=gcbf;
    end
else
    curfig=varargin{1};
end
%%
