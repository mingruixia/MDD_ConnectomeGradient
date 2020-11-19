function varargout = SeeCAT_Montage(varargin)
% SEECAT_MONTAGE MATLAB code for SeeCAT_Montage.fig
%      SEECAT_MONTAGE, by itself, creates a new SEECAT_MONTAGE or raises the existing
%      singleton*.
%
%      H = SEECAT_MONTAGE returns the handle to a new SEECAT_MONTAGE or the handle to
%      the existing singleton*.
%
%      SEECAT_MONTAGE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEECAT_MONTAGE.M with the given input arguments.
%
%      SEECAT_MONTAGE('Property','Value',...) creates a new SEECAT_MONTAGE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SeeCAT_Montage_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SeeCAT_Montage_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SeeCAT_Montage

% Last Modified by GUIDE v2.5 20-Jun-2016 12:46:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SeeCAT_Montage_OpeningFcn, ...
                   'gui_OutputFcn',  @SeeCAT_Montage_OutputFcn, ...
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


% --- Executes just before SeeCAT_Montage is made visible.
function SeeCAT_Montage_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SeeCAT_Montage (see VARARGIN)
if nargin<2+3 && ishandle(varargin{1})
    curfig=w_Compatible2014bFig(varargin{1});
    global st
    handles.st=st{curfig};
    handles.MainFig=curfig;
else
    handles.MainFig=[];
end
%% Init Para
MontageType='C';
Index=-23:5:12;
M=1;
N=8;
I=Inf(M*N, 1);
I(1:length(Index))=Index;
I=reshape(I, [N, M]);
I=I';

%% Init GUI
handles.ImageView=image('Parent', handles.MontageAxe);
handles.lx = line(0, 0, 'Parent',handles.MontageAxe, 'Visible', 'Off');
handles.ly = line(0, 0, 'Parent',handles.MontageAxe, 'Visible', 'Off');

set(handles.MontageAxe, 'XTick', []);
set(handles.MontageAxe, 'XTickLabel', []);
set(handles.MontageAxe, 'YTick', []);
set(handles.MontageAxe, 'YTickLabel', []);
set(handles.MontageAxe, 'YDir', 'Normal');

colormap(handles.MontageAxe, 'gray(64)');
handles=MontageImage(I, MontageType, handles);
%UpdatePosEntry(handles);
RedrawXhairs(handles);

set(handles.ImageView, 'HitTest', 'off', 'Cdata', handles.Image);
set(handles.MontageAxe, 'XLim', [1, size(handles.Image,2)]);
set(handles.MontageAxe, 'YLim', [1, size(handles.Image,1)]);
axis(handles.MontageAxe, 'image');

ResizeFig(handles, size(handles.Image, 2)/size(handles.Image, 1));

set(handles.REntry, 'String', num2str(M));
set(handles.LEntry, 'String', num2str(N));
set(handles.BeginEntry, 'String', num2str(I(1)));
set(handles.DeltasEntry, 'String', '5');
% Choose default command line output for SeeCAT_Montage
handles.IndexM=I;
handles.MontageType=MontageType;
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SeeCAT_Montage wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SeeCAT_Montage_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function ResizeFig(handles, scale)
hold on
set(handles.figure1, 'Unit', 'Pixels');
set(handles.MainFrame, 'Unit', 'Pixels');
set(handles.MontageAxe, 'Unit', 'Pixels');

% Main Frame width height pixel;
posMain=get(handles.MainFrame, 'Position');
xMain=posMain(1);
yMain=posMain(2);
widthMain=posMain(3);
heightMain=posMain(4);
% Gap
heightFull=heightMain./0.3;
gapTop=heightFull*0.025;
gapMiddle=heightFull*0.031;
gapBottom=heightFull*0.02;

widthMontage=widthMain;
heightMontage=widthMontage/scale;
% Fig Pos
heightFig=gapBottom+heightMontage+gapMiddle+heightMain+gapTop;
widthFig=widthMain/0.99;
posFig=get(handles.figure1, 'Position');
posFig(3)=widthFig;
posFig(4)=heightFig;
set(handles.figure1, 'Position', posFig);

% Main Pos
newposMain=posMain;
newposMain(1)=xMain;
newposMain(2)=gapBottom+heightMontage+gapMiddle;
set(handles.MainFrame, 'Position', newposMain);

% Montage Pos
newposMontage(1)=xMain;
newposMontage(2)=gapBottom;
newposMontage(3)=widthMontage;
newposMontage(4)=heightMontage;
set(handles.MontageAxe, 'Position', newposMontage);

%Normalized
set(handles.figure1, 'Unit', 'Normalized');
set(handles.MainFrame, 'Unit', 'Normalized');
set(handles.MontageAxe, 'Unit', 'Normalized');

hold off

function REntry_Callback(hObject, eventdata, handles)
% hObject    handle to REntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of REntry as text
%        str2double(get(hObject,'String')) returns contents of REntry as a double


% --- Executes during object creation, after setting all properties.
function REntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to REntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LEntry_Callback(hObject, eventdata, handles)
% hObject    handle to LEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LEntry as text
%        str2double(get(hObject,'String')) returns contents of LEntry as a double


% --- Executes during object creation, after setting all properties.
function LEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BeginEntry_Callback(hObject, eventdata, handles)
% hObject    handle to BeginEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BeginEntry as text
%        str2double(get(hObject,'String')) returns contents of BeginEntry as a double


% --- Executes during object creation, after setting all properties.
function BeginEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BeginEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DeltasEntry_Callback(hObject, eventdata, handles)
% hObject    handle to DeltasEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DeltasEntry as text
%        str2double(get(hObject,'String')) returns contents of DeltasEntry as a double


% --- Executes during object creation, after setting all properties.
function DeltasEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DeltasEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in IndexButton.
function IndexButton_Callback(hObject, eventdata, handles)
% hObject    handle to IndexButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
options.Resize='on';
options.WindowStyle='modal';

I=inputdlg('Index Matrix', 'Set Index Matrix',...
    size(handles.IndexM, 1)*3,...
    {num2str(handles.IndexM)}, options);
if isempty(I)
    return
end

I=str2num(I{1});
[Image, Space]=w_MontageImage(I, handles.MontageType, handles.MainFig);
set(handles.ImageView, 'HitTest', 'off', 'Cdata', Image);
set(handles.MontageAxe, 'XLim', [1, size(Image,2)]);
set(handles.MontageAxe, 'YLim', [1, size(Image,1)]);

ResizeFig(handles, size(Image, 2)/size(Image, 1));

handles.IndexM=I;
handles.Space=Space;
guidata(hObject, handles);

% --- Executes on button press in ApplyButton.
function ApplyButton_Callback(hObject, eventdata, handles)
% hObject    handle to ApplyButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=str2double(get(handles.REntry, 'String'));
N=str2double(get(handles.LEntry, 'String'));
Begin=str2double(get(handles.BeginEntry, 'String'));
deltas=str2double(get(handles.DeltasEntry, 'String'));
End=Begin+M*N*deltas;
if deltas>0
    End=End-1;
else
    End=End+1;
end

Index=Begin:deltas:End;
I=Inf(M*N, 1);
I(1:end)=Index;
I=reshape(I, [N, M]);
I=I';

if ~isempty(handles.MainFig)
    if ishandle(handles.MainFig)
        global st
        curfig=w_Compatible2014bFig(handles.MainFig);
        handles.st=st{curfig};
    else
        handles.MainFig=[];
    end
end

handles=MontageImage(I, handles.MontageType, handles);
set(handles.ImageView, 'HitTest', 'off', 'Cdata', handles.Image);
set(handles.MontageAxe, 'XLim', [1, size(handles.Image,2)]);
set(handles.MontageAxe, 'YLim', [1, size(handles.Image,1)]);

ResizeFig(handles, size(handles.Image, 2)/size(handles.Image, 1));

handles.IndexM=I;

guidata(hObject, handles);
RedrawXhairs(handles);

% --- Executes on button press in SaveButton.
function SaveButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[File, Path] = uiputfile({'*.tiff';'*.jpeg';'*.png';'*.bmp'},...
    'Save Picture As');
if ~ischar(File)
    return;
end
[Path, Name, Ext]=fileparts(fullfile(Path, File));
IName=sprintf('%s_Index', Name);

%Data=getframe(handles.MontageAxe);
Data=flipdim(get(handles.ImageView, 'CData'), 1);
Data=imresize(Data, 5*[size(Data, 1), size(Data, 2)]);
imwrite(Data, fullfile(Path, [Name, Ext]));
eval(['print -r300 -dtiff -noui ''',fullfile(Path, [Name, '_300dpi', Ext]),''';']); %YAN Chao-Gan, 140806.

if isunix
    Line='unix';
else
    Line='pc';
end
dlmwrite(fullfile(Path, [IName, '.txt']),...
    handles.IndexM, 'delimiter', '\t', 'newline', Line, 'precision', '%d');


% --- Executes on button press in CrosshairButton.
function CrosshairButton_Callback(hObject, eventdata, handles)
% hObject    handle to CrosshairButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Value=get(handles.CrosshairButton, 'Value');
if Value
    State='On';
else
    State='Off';
end
RedrawXhairs(handles);
set(handles.lx, 'Visible', State);
set(handles.ly, 'Visible', State);
% Hint: get(hObject,'Value') returns toggle state of CrosshairButton

function RedrawXhairs(handles)
if ~get(handles.CrosshairButton, 'Value')
    return;
end

st=handles.st;
if nargin<3
    pos=st.centre(:);
end    
bb=st.bb;

Dims = round(diff(bb)'+1);
is   = inv(handles.Space);
cent = is(1:3,1:3)*pos + is(1:3,4);

Type=handles.MontageType;

IndexM=handles.IndexM;
IndexM=flipdim(IndexM, 1);

switch upper(Type)
    case 'T'
        Slice=round(pos(3));
        offsetX=cent(2)-bb(1,2)+1;
        offsetY=cent(1)-bb(1,1)+1;
        D = Dims([1 2]);
    case 'C'
        Slice=round(pos(2));
        offsetX=cent(3)-bb(1,3)+1;
        offsetY=cent(1)-bb(1,1)+1;
        D = Dims([1 3]);
    case 'S'
        Slice=round(pos(1));
        if st.mode==0
            offsetX=cent(2)-bb(1,2)+1;
            offsetY=cent(3)-bb(1,3)+1;
            D = Dims([3 2]);
        else
            offsetX=cent(3)-bb(1,3)+1;
            offsetY=bb(2,2)+1-cent(2);
            D = Dims([2 3]);
        end
end

index=find(IndexM==Slice);
if isempty(index)
    set(handles.lx, 'Xdata', [], 'Ydata', []);
    set(handles.ly, 'Xdata', [], 'Ydata', []);
    return
end

[y, x]=ind2sub(size(IndexM), index);
set(handles.lx,'HitTest','off',...
    'Xdata', [(x-1)*D(1) x*D(1)-1]+0.5,...
    'Ydata', [1 1]*(offsetX+(y-1)*D(2)));
set(handles.ly,'HitTest','off',...
    'Ydata', [(y-1)*D(2) y*D(2)-1]+0.5,...
    'Xdata', [1 1]*(offsetY+(x-1)*D(1)));

% --- Executes on mouse press over axes background.
function MontageAxe_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to MontageAxe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fig=handles.figure1;
if ~strcmpi(get(fig,'SelectionType'),'alt')
    set(fig,...
        'windowbuttonmotionfcn', @(objH, eventData) MontageReposMove(objH, eventData, fig),...
        'windowbuttonupfcn',@(objH, eventData) MontageReposEnd(objH, eventData, fig));
    MontageRepos(fig);
end

function MontageReposMove(objH, eventData, fig)
MontageRepos(fig);

function MontageReposEnd(objH, eventData, fig)
set(fig,'windowbuttonmotionfcn','', 'windowbuttonupfcn','');

function MontageRepos(fig)
handles=guidata(fig);
st=handles.st;
bb=st.bb;
Dims = round(diff(bb)'+1);

handles=guidata(fig);
IndexM=handles.IndexM;
IndexM=flipdim(IndexM, 1);
cp=get(handles.MontageAxe, 'CurrentPoint');
cp=cp(1,1:2);

Type=handles.MontageType;
switch upper(Type)
    case 'T'
        D=Dims([1 2]);
    case 'C'
        D=Dims([1 3]);
    case 'S'
        if st.mode==0
            D=Dims([3 2]);
        else
            D=Dims([2 3]); 
        end 
end

N=fix(cp(1)/D(1));
M=fix(cp(2)/D(2));
Slice=IndexM(M+1, N+1);

cp(1)=cp(1)-N*D(1);
cp(2)=cp(2)-M*D(2);

Space=handles.Space;
is   = inv(Space);
cent = is(1:3,1:3)*zeros(3,1) + is(1:3,4);
switch upper(Type)
    case 'T'
        cent([1 2])=[cp(1)+bb(1,1)-1 cp(2)+bb(1,2)-1];
        centre=Space(1:3,1:3)*cent(:) +Space(1:3,4);
        centre(3)=Slice;
    case 'C'
        cent([1 3])=[cp(1)+bb(1,1)-1 cp(2)+bb(1,3)-1];
        centre=Space(1:3,1:3)*cent(:) +Space(1:3,4);
        centre(2)=Slice;
    case 'S'
        if st.mode ==0
        	cent([3 2])=[cp(1)+bb(1,3)-1 cp(2)+bb(1,2)-1];
        else
            cent([2 3])=[bb(2,2)+1-cp(1) cp(2)+bb(1,3)-1];
        end
        centre=Space(1:3,1:3)*cent(:) +Space(1:3,4);
        centre(1)=Slice;
end

%y_spm_orthviews('Reposition', centre);
handles.st.centre=centre;
guidata(fig, handles)
RedrawXhairs(handles);
if ~isempty(handles.MainFig);
    if ishandle(handles.MainFig)
        curfig=w_Compatible2014bFig(handles.MainFig);
        w_spm_orthviews('Reposition', centre, curfig);
    else
        handles.MainFig=[];
        guidata(fig, handles);
    end
end

% --- Executes on button press in ReduceButton.
function ReduceButton_Callback(hObject, eventdata, handles)
% hObject    handle to ReduceButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Increment=str2double(get(handles.DeltasEntry, 'String'));
if isnan(Increment)
    return
end
handles.IndexM=handles.IndexM-Increment;
set(handles.BeginEntry, 'String', num2str(handles.IndexM(1,1)));

handles=MontageImage(handles.IndexM, handles.MontageType, handles);
set(handles.ImageView, 'HitTest', 'off', 'Cdata', handles.Image);
set(handles.MontageAxe, 'XLim', [1, size(handles.Image,2)]);
set(handles.MontageAxe, 'YLim', [1, size(handles.Image,1)]);

%ResizeFig(handles, size(handles.Image, 2)/size(handles.Image, 1));

guidata(hObject, handles);
RedrawXhairs(handles);

% --- Executes on button press in PlusButton.
function PlusButton_Callback(hObject, eventdata, handles)
% hObject    handle to PlusButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Increment=str2double(get(handles.DeltasEntry, 'String'));
if isnan(Increment)
    return
end
handles.IndexM=handles.IndexM+Increment;
set(handles.BeginEntry, 'String', num2str(handles.IndexM(1,1)));

handles=MontageImage(handles.IndexM, handles.MontageType, handles);
set(handles.ImageView, 'HitTest', 'off', 'Cdata', handles.Image);
set(handles.MontageAxe, 'XLim', [1, size(handles.Image,2)]);
set(handles.MontageAxe, 'YLim', [1, size(handles.Image,1)]);

%ResizeFig(handles, size(handles.Image, 2)/size(handles.Image, 1));

guidata(hObject, handles);
RedrawXhairs(handles);

function handles=AddOverlay(handles, vol)
%%
st=handles.st;
mat = vol.mat;
colourmap=vol.ColorMap;
prop=1-vol.Transparency;
mx=vol.PMax;
mn=vol.NMax;

if ~isfield(st.vols{1},'blobs')
    st.vols{1}.blobs=cell(1,1);
    bset = 1;
else
    bset = numel(st.vols{1}.blobs)+1;
end

c = struct('cmap', colourmap,'prop',prop);
st.vols{1}.blobs{bset} = struct('vol',vol, 'mat',mat, ...
    'max',mx, 'min',mn, 'colour',c);
%addcolourbar(1,bset);
st.curblob=bset;

%Add by Sandy to change colourbar when add overlay
if ~isfield(st.vols{1}.blobs{bset},'colour')
    cmap = get(st.fig,'Colormap');
    if size(cmap,1)~=128
        figure(st.fig)
        spm_figure('Colormap','gray-hot')
    end
    %redraw_colourbar(1,bset,[mn mx],(1:64)'+64);
elseif isstruct(st.vols{1}.blobs{bset}.colour)
    csz   = size(st.vols{1}.blobs{bset}.colour.cmap);
    cdata = reshape(st.vols{1}.blobs{bset}.colour.cmap, [csz(1) 1 csz(2)]);
    %redraw_colourbar(1,bset,[mn mx],cdata);
end
handles.st=st;

function handles=MontageImage(IndexM, type, handles) 
% Written by Xindi WANG
% State Key Laboratory of Cognitive Neuroscience and Learning & IDG/McGovern 
% Institute for Brain Research, Beijing Normal University, Beijing, China
% sandywang.rest@gmail.com
%==========================================================================
st=handles.st;

IndexM=flipdim(IndexM, 1);

bb=st.bb;
Dims = round(diff(bb)'+1);
Space = st.Space;

switch upper(type)
    case 'T'
        R=Dims(2)*size(IndexM, 1);
        L=Dims(1)*size(IndexM, 2);
    case 'C'
        R=Dims(3)*size(IndexM, 1);
        L=Dims(1)*size(IndexM, 2);
    case 'S'
        if st.mode == 0
            R=Dims(2)*size(IndexM, 1);
            L=Dims(3)*size(IndexM, 2);
        else
            R=Dims(3)*size(IndexM, 1);
            L=Dims(2)*size(IndexM, 2);
        end
end
if isfield(st.vols{1},'blobs')
    Image=zeros(R, L, 3);
else
    Image=zeros(R, L);
end

for index=1:numel(IndexM)
    [y, x]=ind2sub(size(IndexM), index);
    slice=IndexM(y, x);
    if isinf(slice)
        continue;
    end
    is=inv(st.Space);
    
    switch upper(type)
        case 'T'
            cent = is(1:3,1:3)*[0;0;slice] + is(1:3,4);
            M = st.Space\st.vols{1}.premul*st.vols{1}.mat;
            AM0 = [ 1 0 0 -bb(1,1)+1
                    0 1 0 -bb(1,2)+1
                    0 0 1 -cent(3)
                    0 0 0 1];
            AM = inv(AM0*M);
            AD = Dims([1 2]);
        case 'C'
            cent = is(1:3,1:3)*[0;slice;0] + is(1:3,4);
            M = st.Space\st.vols{1}.premul*st.vols{1}.mat;
            AM0 = [ 1 0 0 -bb(1,1)+1
                    0 0 1 -bb(1,3)+1
                    0 1 0 -cent(2)
                    0 0 0 1];
            AM = inv(AM0*M);
            AD = Dims([1 3]);
        case 'S'
            cent = is(1:3,1:3)*[slice;0;0] + is(1:3,4);
            M = st.Space\st.vols{1}.premul*st.vols{1}.mat;
            if st.mode ==0
                AM0 = [ 0 0 1 -bb(1,3)+1
                        0 1 0 -bb(1,2)+1
                        1 0 0 -cent(1)
                        0 0 0 1];
                AM = inv(AM0*M); 
                AD = Dims([3 2]);
            else
                AM0 = [ 0 -1 0 +bb(2,2)+1
                        0  0 1 -bb(1,3)+1
                        1  0 0 -slice
                        0  0 0 1];
                AM = inv(AM0*M);
                AD = Dims([2 3]);
            end
    end
    
    try
        %Revised by YAN Chao-Gan, 130720. Could also work with Data has been read into memory other than only depending on the file.
        if ~isfield(st.vols{1},'Data')
            imga = spm_slice_vol(st.vols{1},AM,AD,st.hld)';
        else   %Revised by YAN Chao-Gan, 130720. Could also work with Data has been read into memory other than only depending on the file.
            imga = spm_slice_vol(st.vols{1}.Data,AM,AD,st.hld)';
        end
        
        ok   = true;
    catch
        fprintf('Cannot access file "%s".\n', st.vols{1}.fname);
        fprintf('%s\n',getfield(lasterror,'message'));
        ok   = false;
    end
    if ok
        % get min/max threshold
        if strcmp(st.vols{1}.window,'auto')
            mn = -Inf;
            mx = Inf;
        else
            mn = min(st.vols{1}.window);
            mx = max(st.vols{1}.window);
        end
        % threshold images
        imga = max(imga,mn); imgt = min(imga,mx);
        % compute intensity mapping, if histeq is available
        if license('test','image_toolbox') == 0
            st.vols{1}.mapping = 'linear';
        end
        switch st.vols{1}.mapping
            case 'linear'
            case 'histeq'
                % scale images to a range between 0 and 1
                imga1=(imga-min(imga(:)))/(max(imga(:)-min(imga(:)))+eps);
                img  = histeq(imga1(:),1024);
                imga = reshape(img(1:numel(imga1)),size(imga1));
                mn = 0;
                mx = 1;
            case 'quadhisteq'
                % scale images to a range between 0 and 1
                imga1=(imga-min(imga(:)))/(max(imga(:)-min(imga(:)))+eps);
                img  = histeq(imga1(:).^2,1024);
                imga = reshape(img(1:numel(imga1)),size(imga1));
                mn = 0;
                mx = 1;
            case 'loghisteq'
                sw = warning('off','MATLAB:log:logOfZero');
                imga = log(imga-min(imga(:)));
                warning(sw);
                imga(~isfinite(imga)) = 0;
                % scale log images to a range between 0 and 1
                imga1=(imga-min(imga(:)))/(max(imga(:)-min(imga(:)))+eps);
                img  = histeq(imga1(:),1024);
                imga = reshape(img(1:numel(imga1)),size(imga1));
                mn = 0;
                mx = 1;
        end
        % recompute min/max for display
        if strcmp(st.vols{1}.window,'auto')
            mx = -inf; mn = inf;
        end
        %Add by Sandy, make the same mn/mx in a volume
        if ~isfield(st.vols{1},'Data')
            if ~isempty(imga)
                tmp = imgt(isfinite(imga));
                mx = max([mx max(max(tmp))]);
                mn = min([mn min(min(tmp))]);
            end
        else
            mx=max(max(max(st.vols{1}.Data)));
            mn=min(min(min(st.vols{1}.Data)));
        end
        if mx==mn, mx=mn+eps; end
        
        if isfield(st.vols{1},'blobs')
            if isstruct(st.vols{1}.blobs{1}.colour)
                % Add blobs for display using a defined colourmap
                
                % colourmaps
                gryc = (0:63)'*ones(1,3)/63;
                
                % scale grayscale image, not isfinite -> black
                imga = scaletocmap(imga,mn,mx,gryc,65);

                gryc = [gryc; 0 0 0];
                
                % combine gray and blob data to
                % truecolour
                
                %Modified by Sandy for Multi-Overlay 20140104
                umpa=ones(size(imga));
                umpa=repmat(umpa(:),1,3);
                ompa=zeros(size(umpa));
                for j=1:numel(st.vols{1}.blobs)
                    actc = st.vols{1}.blobs{j}.colour.cmap;
                    actp = st.vols{1}.blobs{j}.colour.prop;
                    % get max for blob image
                    if isfield(st.vols{1}.blobs{j},'max')
                        cmx = st.vols{1}.blobs{j}.max;
                    else
                        cmx = max([eps maxval(st.vols{1}.blobs{j}.vol)]);
                    end
                    if isfield(st.vols{1}.blobs{j},'min')
                        cmn = st.vols{1}.blobs{j}.min;
                    else
                        cmn = -cmx;
                    end
                    
                    % get blob data
                    vol  = st.vols{1}.blobs{j}.vol;
                    M    = st.Space\st.vols{1}.premul*st.vols{1}.blobs{j}.mat;
                
                    %Revised by YAN Chao-Gan, 130720. Could also work with Data has been read into memory other than only depending on the file.
                    if ~isfield(vol,'Data')
                        tmpa = spm_slice_vol(vol,inv(AM0*M),AD,[0 NaN])';
                    else   %Revised by YAN Chao-Gan, 130720. Could also work with Data has been read into memory other than only depending on the file.
                        tmpa = spm_slice_vol(vol.Data,inv(AM0*M),AD,[0 NaN])';
                    end
                
                
                    % actimg scaled round 0, black NaNs
                    topc = size(actc,1)+1;
                    tmpa = scaletocmap(tmpa,cmn,cmx,actc,topc);
              
                    %Overlay Transparent Weight
                    jmpa = actp*repmat((tmpa(:)~=topc),1,3);
                    
                    %Except Underlay
                    umpa = umpa-jmpa;
                    
                    %Negtive Weight Recoup
                    nmpa = umpa.*(umpa<0);
                    
                    %Modified Overlay Transparent Weight
                    jmpa=jmpa+nmpa;
                    
                    actc = [actc; 0 0 0];
                
                    ompa = ompa+jmpa.*actc(tmpa(:),:);
                    
                    umpa(umpa<0) = 0;
                end      
                imga = reshape(ompa+gryc(imga(:),:).*umpa, [size(imga) 3]);
                Image(AD(2)*(y-1)+1:AD(2)*y, AD(1)*(x-1)+1:AD(1)*x, :)=imga;
            end
        else
            scal = 64/(mx-mn);
            dcoff = -mn*scal;
            imga = imga*scal+dcoff;
            Image(AD(2)*(y-1)+1:AD(2)*y, AD(1)*(x-1)+1:AD(1)*x)=imga;
        end        
    end
end

handles.Image=Image;
handles.Space=Space;

function img = scaletocmap(inpimg,mn,mx,cmap,miscol)
if nargin < 5, miscol=1; end
cml = size(cmap,1);
scf = (cml-1)/(mx-mn);
img = round((inpimg-mn)*scf)+1;
img(img<1)   = 1; 
img(img>cml) = cml;
img(inpimg==0) = miscol; %Added by YAN Chao-Gan 130609, mask out the 0 voxels.
img(~isfinite(img)) = miscol;



% --- Executes on selection change in MontageTypePopup.
function MontageTypePopup_Callback(hObject, eventdata, handles)
% hObject    handle to MontageTypePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Value=get(handles.MontageTypePopup, 'Value');
switch upper(Value)
    case 1 %C
        MontageType='C';
    case 2 %S
        MontageType='S';
    case 3 %T
        MontageType='T';
end
I=handles.IndexM;
handles=MontageImage(I, MontageType, handles);
%UpdatePosEntry(handles);
handles.MontageType=MontageType;
RedrawXhairs(handles);
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns MontageTypePopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MontageTypePopup


% --- Executes during object creation, after setting all properties.
function MontageTypePopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MontageTypePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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


% --- Executes on button press in CIReportBtn.
function CIReportBtn_Callback(hObject, eventdata, handles)
% hObject    handle to CIReportBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
