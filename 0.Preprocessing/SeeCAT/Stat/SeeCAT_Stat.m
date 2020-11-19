function varargout = SeeCAT_Stat(varargin)
% SeeCAT_Stat MATLAB code for SeeCAT_Stat.fig
%      SeeCAT_Stat, by itself, creates a new SeeCAT_Stat or raises the existing
%      singleton*.
%
%      H = SeeCAT_Stat returns the handle to a new SeeCAT_Stat or the handle to
%      the existing singleton*.
%
%      SeeCAT_Stat('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SeeCAT_Stat.M with the given input arguments.
%
%      SeeCAT_Stat('Property','Value',...) creates a new SeeCAT_Stat or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SeeCAT_Stat_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SeeCAT_Stat_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SeeCAT_Stat

% Last Modified by GUIDE v2.5 22-Jun-2016 21:27:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SeeCAT_Stat_OpeningFcn, ...
                   'gui_OutputFcn',  @SeeCAT_Stat_OutputFcn, ...
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
%% Revided from DPABI_STAT_TOOL by Sandy Wang

% --- Executes just before SeeCAT_Stat is made visible.
function SeeCAT_Stat_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SeeCAT_Stat (see VARARGIN)
if nargin > 3
    Flag=varargin{1};
else
    Flag='T1';
end
[handles.SampleNum, Value]=SwitchType(Flag);
set(handles.TypeListbox, 'Value', Value);
handles=StatType(handles);
%handles=ClearConfigure(handles);
handles.SampleCells={};
handles.CovImageCells={};
handles.CurDir=pwd;
set(handles.OutputDirEntry, 'String', pwd);
% Choose default command line output for SeeCAT_Stat
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SeeCAT_Stat wait for user response (see UIRESUME)
% uiwait(handles.figure1);
ProgramPath = fileparts(which('SeeCAT.m'));
if ispc
    filesep = '\';
else
    filesep = '/';
end
set(handles.MaskEntry, 'String',[ProgramPath,filesep,'templates',filesep,'GreyMask_02_61x73x61.nii']);

function [Lim, Value]=SwitchType(Flag)
switch upper(Flag)
    case 'T1'
        Value=1;     
        Lim=1;
    case 'T2'
        Value=2;
        Lim=2;
    case 'TP'
        Value=3;
        Lim=2;
    case 'F'
        Value=4;
        Lim=100;
    case 'FR'
        Value=5;
        Lim=100;
    case 'R'
        Value=6;
        Lim=1;
end

function handles=StatType(handles)
Value=get(handles.TypeListbox, 'Value');
CorrFlag='Off';
BaseFlag='Off';

switch Value
    case 1
        handles.SampleNum=1;
        Prefix='T1';
        BaseFlag='On';
    case 2
        handles.SampleNum=2;
        Prefix='T2';
    case 3
        handles.SampleNum=2;
        Prefix='TP';
    case 4
        handles.SampleNum=10;
        Prefix='F';
    case 5
        handles.SampleNum=10;
        Prefix='FR';
    case 6
        handles.SampleNum=1;
        Prefix='R';
        CorrFlag='On';
end

handles=ClearConfigure(handles);

set(handles.BaseLab,   'Visible', BaseFlag);
set(handles.BaseEntry, 'Visible', BaseFlag);
set(handles.PrefixEntry, 'String', Prefix);
%set(handles.CorrSeedListbox,      'Visible', CorrFlag);
%set(handles.CorrSeedRemoveButton, 'Visible', CorrFlag);
%set(handles.CorrSeedAddButton,    'Enable', CorrFlag);
set(handles.CorrSeedListbox,        'Enable', CorrFlag); 
set(handles.CorrSeedAddButton,      'Enable', CorrFlag); 
set(handles.CorrSeedRemoveButton,   'Enable', CorrFlag); 


% --- Outputs from this function are returned to the command line.
function varargout = SeeCAT_Stat_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in CovTextListbox.
function CovTextListbox_Callback(hObject, eventdata, handles)
% hObject    handle to CovTextListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CovTextListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CovTextListbox


% --- Executes during object creation, after setting all properties.
function CovTextListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CovTextListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in CovImageListbox.
function CovImageListbox_Callback(hObject, eventdata, handles)
% hObject    handle to CovImageListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CovImageListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CovImageListbox


% --- Executes during object creation, after setting all properties.
function CovImageListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CovImageListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in SampleListbox.
function SampleListbox_Callback(hObject, eventdata, handles)
% hObject    handle to SampleListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SampleListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SampleListbox


% --- Executes during object creation, after setting all properties.
function SampleListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SampleListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CovTextAddButton.
function CovTextAddButton_Callback(hObject, eventdata, handles)
% hObject    handle to CovTextAddButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Name, Path]=uigetfile({'*.txt;*.csv;*.tsv','Text Covariates File (*.txt;*.csv;*.tsv)';'*.*', 'All Files (*.*)';},...
    'Pick the Text Covariates', 'MultiSelect','on', handles.CurDir);
if isnumeric(Name)
    return
end

if ischar(Name)
    Name={Name};
end
Name=Name';
PathCell=cellfun(@(name) fullfile(Path, name), Name, 'UniformOutput', false);
AddString(handles.CovTextListbox, PathCell);
UpdatePreview(handles);

% --- Executes on button press in CovTextRemoveButton.
function CovTextRemoveButton_Callback(hObject, eventdata, handles)
% hObject    handle to CovTextRemoveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Value=get(handles.CovTextListbox, 'Value');
if Value==0
    return
end
RemoveString(handles.CovTextListbox, Value);
UpdatePreview(handles);

% --- Executes on button press in CovImageAddButton.
function CovImageAddButton_Callback(hObject, eventdata, handles)
% hObject    handle to CovImageAddButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Path=uigetdir(handles.CurDir, 'Pick Covariate Images Directory');
if isnumeric(Path)
    return
end
[handles.CurDir, Name]=fileparts(Path);

D=dir(fullfile(Path, '*.hdr'));
if isempty(D)
    D=dir(fullfile(Path, '*.nii'));
end

if isempty(D)
    D=dir(fullfile(Path, '*.nii.gz'));
end
NameCell={D.name}';
Num=numel(NameCell);
CovImageCell=cellfun(@(Name) fullfile(Path, Name), NameCell,...
    'UniformOutput', false);
handles.CovImageCells{numel(handles.CovImageCells)+1}=CovImageCell;
StringOne={sprintf('[%d] (%s) %s', Num, Name, Path)};
AddString(handles.CovImageListbox, StringOne);
guidata(hObject, handles);
UpdatePreview(handles);

% --- Executes on button press in CovImageRemoveButton.
function CovImageRemoveButton_Callback(hObject, eventdata, handles)
% hObject    handle to CovImageRemoveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Value=get(handles.CovImageListbox, 'Value');
if Value==0
    return
end
handles.CovImageCells(Value)=[];
RemoveString(handles.CovImageListbox, Value);
guidata(hObject, handles);
UpdatePreview(handles);

% --- Executes on button press in SampleAddButton.
function SampleAddButton_Callback(hObject, eventdata, handles)
% hObject    handle to SampleAddButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if numel(handles.SampleCells)==handles.SampleNum
    warndlg('Invalid Number of Groups');
    return
end
Path=uigetdir(handles.CurDir, 'Pick Sample Directory');
if isnumeric(Path)
    return
end
[handles.CurDir, Name]=fileparts(Path);

D=dir(fullfile(Path, '*.hdr'));
if isempty(D)
    D=dir(fullfile(Path, '*.nii'));
end

if isempty(D)
    D=dir(fullfile(Path, '*.nii.gz'));
end
NameCell={D.name}';
Num=numel(NameCell);
SampleCell=cellfun(@(Name) fullfile(Path, Name), NameCell,...
    'UniformOutput', false);
handles.SampleCells{numel(handles.SampleCells)+1}=SampleCell;
StringOne={sprintf('[%d] (%s) %s', Num, Name, Path)};
AddString(handles.SampleListbox, StringOne);
guidata(hObject, handles);
UpdatePreview(handles);

% --- Executes on button press in SampleRemoveButton.
function SampleRemoveButton_Callback(hObject, eventdata, handles)
% hObject    handle to SampleRemoveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Value=get(handles.SampleListbox, 'Value');
if Value==0
    return
end
handles.SampleCells(Value)=[];
RemoveString(handles.SampleListbox, Value);
guidata(hObject, handles);
UpdatePreview(handles);

function MaskEntry_Callback(hObject, eventdata, handles)
% hObject    handle to MaskEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaskEntry as text
%        str2double(get(hObject,'String')) returns contents of MaskEntry as a double


% --- Executes during object creation, after setting all properties.
function MaskEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaskEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MaskButton.
function MaskButton_Callback(hObject, eventdata, handles)
% hObject    handle to MaskButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Name, Path]=uigetfile({'*.img;*.nii;*.nii.gz','Brain Image Files (*.img;*.nii;*.nii.gz)';'*.*', 'All Files (*.*)';},...
    'Pick the Mask file', handles.CurDir);
if isnumeric(Name)
    return
end
set(handles.MaskEntry, 'String', fullfile(Path, Name));


function OutputDirEntry_Callback(hObject, eventdata, handles)
% hObject    handle to OutputDirEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OutputDirEntry as text
%        str2double(get(hObject,'String')) returns contents of OutputDirEntry as a double


% --- Executes during object creation, after setting all properties.
function OutputDirEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OutputDirEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in OutputDirButton.
function OutputDirButton_Callback(hObject, eventdata, handles)
% hObject    handle to OutputDirButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Path=uigetdir(handles.CurDir, 'Pick Output Directory');
if isnumeric(Path)
    return
end
handles.CurDir=fileparts(Path);
set(handles.OutputDirEntry, 'String', Path);
guidata(hObject, handles);


function PrefixEntry_Callback(hObject, eventdata, handles)
% hObject    handle to PrefixEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PrefixEntry as text
%        str2double(get(hObject,'String')) returns contents of PrefixEntry as a double


% --- Executes during object creation, after setting all properties.
function PrefixEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PrefixEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ComputeButton.
function ComputeButton_Callback(hObject, eventdata, handles)
% hObject    handle to ComputeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.SampleCells)
    return
end

[Y_Volume, C_Volume, C_Text, R_Text, NiiHeader, NumOfSamples]=GetVolume(handles);
Vox=sqrt(sum(NiiHeader.mat(1:3,1:3).^2));
SumNum=sum(NumOfSamples);

Value=get(handles.TypeListbox, 'Value');
DF_G=[];

switch Value
    case 1 %T1
        Base=str2double(get(handles.BaseEntry, 'String'));
        if isnan(Base)
            errordlg('Invalid Base Value');
            return
        end
        
        Y_Volume=Y_Volume-Base;
        if ~isempty(C_Volume)
            C_Volume=C_Volume-repmat(mean(C_Volume, 4), [1,1,1,SumNum]);
        end
        if ~isempty(C_Text)
            C_Text=C_Text-repmat(mean(C_Text,1),[SumNum,1]);
        end
        Regressor=ones(SumNum, 1);
        
        % Contrast
        Contrast=zeros(1, size(Regressor, 2)+size(C_Text, 2));
        if ~isempty(C_Volume)
            Contrast=[Contrast, 0];
        end        
        Contrast(1)=1;
        TFR_Flag='T';
    case 2 %T2
        GL_Cell=cell(numel(NumOfSamples), 1);
        for i=1:numel(GL_Cell)
            GL_Cell{i, 1}=i*ones(NumOfSamples(i), 1);
        end
        GroupLabel=cell2mat(GL_Cell);        
        GroupLabel(GroupLabel==2)=-1;
        M_Regressor=GroupLabel;
        Regressor=[M_Regressor, ones(SumNum, 1)];
        
        % Contrast
        Contrast=zeros(1, size(Regressor, 2)+size(C_Text, 2));
        if ~isempty(C_Volume)
            Contrast=[Contrast, 0];
        end          
        Contrast(1)=1;        
        TFR_Flag='T';
    case 3 %TP
        GL_Cell=cell(numel(NumOfSamples), 1);
        for i=1:numel(GL_Cell)
            GL_Cell{i, 1}=i*ones(NumOfSamples(i), 1);
        end
        GroupLabel=cell2mat(GL_Cell);        
        GroupLabel(GroupLabel==2)=-1;
        
        M_Regressor=GroupLabel;
        S_Regressor=zeros(SumNum, NumOfSamples(1));
        for n=1:NumOfSamples(1)
            S_Regressor([n,n + NumOfSamples(1)], n) = 1;
        end
        Regressor=[M_Regressor, S_Regressor];
        
        % Contrast
        Contrast=zeros(1, size(Regressor, 2)+size(C_Text, 2));
        if ~isempty(C_Volume)
            Contrast=[Contrast, 0];
        end          
        Contrast(1)=1;                
        TFR_Flag='T';
    case 4 %F
        GL_Cell=cell(numel(NumOfSamples), 1);
        for i=1:numel(GL_Cell)
            GL_Cell{i, 1}=i*ones(NumOfSamples(i), 1);
        end
        GroupLabel=cell2mat(GL_Cell);        
        
        DF_G=numel(NumOfSamples)-1;
        M_Regressor=zeros(SumNum, DF_G);
        for i=1:DF_G
            M_Regressor(:, i)=GroupLabel==i;
        end
        Regressor=[M_Regressor, ones(SumNum, 1)];
        
        % Contrast
        Contrast=zeros(1, size(Regressor, 2)+size(C_Text, 2));
        if ~isempty(C_Volume)
            Contrast=[Contrast, 0];
        end 
        Contrast(1:DF_G)=1;
        TFR_Flag='F';
    case 5 %FR
        GL_Cell=cell(numel(NumOfSamples), 1);
        for i=1:numel(GL_Cell)
            GL_Cell{i, 1}=i*ones(NumOfSamples(i), 1);
        end
        GroupLabel=cell2mat(GL_Cell);        
        
        DF_G=numel(NumOfSamples)-1;
        M_Regressor=zeros(SumNum, DF_G);
        for i=1:DF_G
            M_Regressor(:, i)=GroupLabel==i;
        end

        S_Regressor=zeros(SumNum, NumOfSamples(1));
        for n=1:NumOfSamples(1)
            S_Regressor(i:n:SumNum, i) = 1;
        end   
        Regressor=[M_Regressor, S_Regressor];
        
        % Contrast
        Contrast=zeros(1, size(Regressor, 2)+size(C_Text, 2));
        if ~isempty(C_Volume)
            Contrast=[Contrast, 0];
        end        
        Contrast(1:DF_G)=1;
        TFR_Flag='F';
    case 6 %R
        RegressorCells=cell(size(R_Text, 2), 1);
        for c=1:size(R_Text, 2)
            M_Regressor=R_Text(:, c);
            Regressor=[M_Regressor, ones(SumNum, 1)];
            RegressorCells{c, 1}=Regressor;
        end
        
        % Contrast
        Contrast=zeros(1, size(Regressor, 2)+size(C_Text, 2));
        if ~isempty(C_Volume)
            Contrast=[Contrast, 0];
        end        
        Contrast(1)=1;                        
        TFR_Flag='R';
end
DF_S=SumNum-numel(Contrast);

% Mask
MaskFile=get(handles.MaskEntry, 'String');
if ~isempty(MaskFile)
    MaskData=y_ReadRPI(MaskFile);
    MaskVector=reshape(MaskData, [], 1);
else
    MaskVector=[];
end

% Data
[N1, N2, N3, N4]=size(Y_Volume);
Y_Volume=reshape(Y_Volume, [], SumNum);
FMaskVector=any(Y_Volume, 2);
if ~isempty(MaskVector)
    FMaskVector=logical(FMaskVector.*MaskVector);   
end
FMaskData=reshape(FMaskVector, [N1, N2, N3]);
Y_Volume=Y_Volume(FMaskVector, :);
if ~isempty(C_Volume)
    C_Volume=reshape(C_Volume, [], SumNum);
    C_Volume=C_Volume(FMaskVector, :);
end

if Value==6 && numel(RegressorCells)>1
    for c=1:numel(RegressorCells)
        [S, R]=Regression(Y_Volume, RegressorCells{c}, C_Text, C_Volume, Contrast, TFR_Flag);
        S_Volume=zeros(N1*N2*N3, 1);
        S_Volume(FMaskVector)=S;
        R_Volume=zeros(N1*N2*N3, N4);
        R_Volume(FMaskVector, :)=R;

        S_Volume=reshape(S_Volume, [N1, N2, N3]);
        R_Volume=reshape(R_Volume, [N1, N2, N3, N4]);
        
        [dLh, FWHM]=Smoothest(R_Volume, FMaskData, DF_S, Vox);
        
        OutputDir=get(handles.OutputDirEntry, 'String');
        if isempty(OutputDir)
            OutputDir=fileparts(handles.CurDir);
        end
        Prefix=get(handles.PrefixEntry, 'String');
        OutputName=fullfile(OutputDir, [Prefix, sprintf('_%d', c),'.nii']);
        ResName=fullfile(OutputDir, [Prefix, sprintf('_%d', c), '_Residual.nii']);
        
        NiiHeader.descrip=sprintf('CAT{T_[%.1f]}{dLh_%f}{FWHMx_%fFWHMy_%fFWHMz_%fmm}',...
            DF_S, dLh, FWHM(1), FWHM(2), FWHM(3));  
        y_Write(S_Volume, NiiHeader, OutputName);
        y_Write(R_Volume, NiiHeader, ResName);        
    end
else
    [S, R]=Regression(Y_Volume, Regressor, C_Text, C_Volume, Contrast, TFR_Flag);
    S_Volume=zeros(N1*N2*N3, 1);
    S_Volume(FMaskVector)=S;
    R_Volume=zeros(N1*N2*N3, N4);
    R_Volume(FMaskVector, :)=R;

    S_Volume=reshape(S_Volume, [N1, N2, N3]);
    R_Volume=reshape(R_Volume, [N1, N2, N3, N4]);

    [dLh, FWHM]=Smoothest(R_Volume, FMaskData, DF_S, Vox);
    
    OutputDir=get(handles.OutputDirEntry, 'String');
    if isempty(OutputDir)
        OutputDir=fileparts(handles.CurDir);
    end
    Prefix=get(handles.PrefixEntry, 'String');
    OutputName=fullfile(OutputDir, [Prefix, '.nii']);
    ResName=fullfile(OutputDir, [Prefix, '_Residual.nii']);
    
    if ~isempty(DF_G)
        NiiHeader.descrip=sprintf('CAT{F_[%.1f,%.1f]}{dLh_%f}{FWHMx_%fFWHMy_%fFWHMz_%fmm}',...
            DF_G,DF_S, dLh, FWHM(1), FWHM(2), FWHM(3));
    else
        NiiHeader.descrip=sprintf('CAT{T_[%.1f]}{dLh_%f}{FWHMx_%fFWHMy_%fFWHMz_%fmm}',...
            DF_S, dLh, FWHM(1), FWHM(2), FWHM(3));        
    end
    y_Write(S_Volume, NiiHeader, OutputName);
    y_Write(R_Volume, NiiHeader, ResName);
end

% --- Executes on button press in HelpButton.
function HelpButton_Callback(hObject, eventdata, handles)
% hObject    handle to HelpButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Value=get(handles.TypeListbox, 'Value');
switch Value
    case 1
        msgbox({'One-Sample T-Test:';...
            'Please specify the group images (Only one group)';...
            'Base means the null hypothesis is that the group mean equals to the base value. Default: 0';...
            'The value of each voxel in the output image is a T statistic value. The degree of freedom information is stored in the header of the output image file.';...
            },'Help');
    case 2
        msgbox({'Two-Sample T-Test:';...
            'If only the group images are specified, then perform voxel-wise Two-Sample T-Test.';...
            'If the covariate images are also specified (e.g. gray matter proportion images), then voxel-wise Two-Sample T-Test is performed while take each voxel in the covariate images as a covaraite. Please make sure the correspondence between the group images and the covariate images.';...
            'Text covariate can be also specified as text files. (E.g. age, brain size, IQ etc.)';...
            'The value of each voxel in the output image is a T statistic value (positive means the mean of Group 1 is greater than the mean of Group 2). The degree of freedom information is stored in the header of the output image file.';...
            },'Help');
    case 3
        msgbox({'Paired T-Test:';...
            'Paired T-Test is performed between the two conditions. Please make sure the correspondence of the images between the two paired conditions.';...
            'The value of each voxel in the output image is a T statistic value (positive means Condition 1 is greater than Condition 2). The degree of freedom information is stored in the header of the output image file.';...
            },'Help');        
    case 4
        msgbox({'ANOVA or ANCOVA analysis:';...
            'If only the group images are specified, then perform voxel-wise ANOVA analysis.';...
            'If the covariate images are also specified (e.g. gray matter proportion images), then voxel-wise ANCOVA analysis is performed while take each voxel in the covariate images as a covaraite. Please make sure the correspondence between the group images and the covariate images.';...
            'Text covariate can be also specified as text files. (E.g. age, brain size, IQ etc.)';...
            'The value of each voxel in the output image is an F statistic value. The degree of freedom information is stored in the header of the output image file.';...
            },'Help');
    case 6
        msgbox({'Correlation Analysis:';...
            'If only the group images and the seed variate are specified, then perform Pearson''s correlation analysis.';...
            'If the covariate images are also specified (e.g. gray matter proportion images), then partial correlation analysis is performed while take each voxel in the covariate images as a covaraite. Please make sure the correspondence between the group images and the covariate images.';...
            'Text covariate can be also specified as text files. (E.g. age, brain size, IQ etc.)';...
            'The value of each voxel in the output image is an R statistic value. The degree of freedom information is stored in the header of the output image file.';...
            },'Help');
end

function handles=ClearConfigure(handles)
set(handles.SampleListbox,   'String', '', 'Value', 0);
handles.SampleCells={};
handles.CovImageCells={};

set(handles.CovImageListbox, 'String', '', 'Value', 0);
set(handles.CovTextListbox,  'String', '', 'Value', 0);
set(handles.CorrSeedListbox, 'String', '', 'Value', 0);

% set(handles.MaskEntry,       'String', '');
set(handles.GLMTable, 'Data', cell(10, 10));
set(handles.GLMTable, 'ColumnName', cell(10, 1));

function AddString(ListboxHandle, NewCell)
StringCell=get(ListboxHandle, 'String');
StringCell=[StringCell; NewCell];
set(ListboxHandle, 'String', StringCell, 'Value', numel(StringCell));

function RemoveString(ListboxHandle, Value)
StringCell=get(ListboxHandle, 'String');
StringCell(Value)=[];
if isempty(StringCell)
    Value=0;
end
if Value > numel(StringCell)
    Value=Value-1;
end
set(ListboxHandle, 'String', StringCell, 'Value', Value);


% --- Executes on selection change in CorrSeedListbox.
function CorrSeedListbox_Callback(hObject, eventdata, handles)
% hObject    handle to CorrSeedListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CorrSeedListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CorrSeedListbox


% --- Executes during object creation, after setting all properties.
function CorrSeedListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CorrSeedListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CorrSeedAddButton.
function CorrSeedAddButton_Callback(hObject, eventdata, handles)
% hObject    handle to CorrSeedAddButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RString=get(handles.CorrSeedAddButton, 'String');
if numel(RString)==handles.SampleNum
    warndlg('Invalid Number of Groups');
    return
end

[Name, Path]=uigetfile({'*.txt;*.csv;*.tsv','Correlation Seed Series File (*.txt;*.csv;*.tsv)';'*.*', 'All Files (*.*)';},...
    'Pick the Text Covariates');
if isnumeric(Name)
    return
end

PathCell={fullfile(Path, Name)};
AddString(handles.CorrSeedListbox, PathCell);

UpdatePreview(handles);

% --- Executes on button press in CorrSeedRemoveButton.
function CorrSeedRemoveButton_Callback(hObject, eventdata, handles)
% hObject    handle to CorrSeedRemoveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Value=get(handles.CorrSeedListbox, 'Value');
if Value==0
    return
end

index=0;
for i=1:numel(handles.SampleCells)
    if ischar(handles.SampleCells{i})
        index=i;
        break
    end
end
RemoveString(handles.CorrSeedListbox, Value);
UpdatePreview(handles);

function BaseEntry_Callback(hObject, eventdata, handles)
% hObject    handle to BaseEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BaseEntry as text
%        str2double(get(hObject,'String')) returns contents of BaseEntry as a double


% --- Executes during object creation, after setting all properties.
function BaseEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BaseEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in TypeListbox.
function TypeListbox_Callback(hObject, eventdata, handles)
% hObject    handle to TypeListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=StatType(handles);

guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns TypeListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TypeListbox


% --- Executes during object creation, after setting all properties.
function TypeListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TypeListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ExportBtn.
function ExportBtn_Callback(hObject, eventdata, handles)
% hObject    handle to ExportBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[File , Path]=uiputfile({'*.txt;*.tsv','Preview Table (*.txt;*.tsv)';'*.*', 'All Files (*.*)';}, ...
    'Export Preview Table' , handles.CurDir);
if ~ischar(File)
    return
end
String=get(handles.GLMTable, 'Data');
ColumnName=get(handles.GLMTable, 'ColumnName');
String=[ColumnName'; String];
fid=fopen(fullfile(Path, File), 'w');
L=size(String, 2);
F=[];
for l=1:L
    F=[F, '%s\t'];
end
F(end)='n';
Tmp=String';
fprintf(fid, F, Tmp{:});
fclose(fid);

function UpdatePreview(handles)
% Cells
SampleCells=handles.SampleCells;
CovImageCells=handles.CovImageCells;
CovTextCells=get(handles.CovTextListbox, 'String');

% Generate String
if ~isempty(SampleCells)
    PreviewString=[];
    ColumnNames=[];
    L=max(cellfun(@numel, SampleCells));
    for i=1:numel(SampleCells)
        ColumnNames=[ColumnNames; {sprintf('Group - %d', i)}];
        SampleString=cellfun(@(f) GetFileName(f), SampleCells{i},...
            'UniformOutput', false);
        TmpString=cell(L, 1);
        TmpString(1:numel(SampleString), 1)=SampleString;
        PreviewString=[PreviewString, TmpString];
        
        %Cov Image
        if ~isempty(CovImageCells) && i<=numel(CovImageCells)
            ColumnNames=[ColumnNames; {sprintf('Cov Image - %d', i)}];
            CovImageString=cellfun(@(f) GetFileName(f), CovImageCells{i},...
                'UniformOutput', false);
            TmpString=cell(L, 1);
            TmpString(1:numel(CovImageString), 1)=CovImageString;            
            PreviewString=[PreviewString, TmpString];
        end
        
        %Cov Text
        if ~isempty(CovTextCells) && i<=numel(CovTextCells)
            CovTextData=load(CovTextCells{i});
            [TextPath, TextName, TextExt]=fileparts(CovTextCells{i});
            NumOfCovText=size(CovTextData, 2);
            for j=1:NumOfCovText
                ColumnNames=[ColumnNames; {sprintf('Cov Text - %s(%d)', TextName, j)}];
            end
            CovTextData=cellfun(@(n) num2str(n), num2cell(CovTextData),...
                'UniformOutput', false);
            TmpString=cell(L, NumOfCovText);
            TmpString(1:size(CovTextData, 1), :)=CovTextData;
            PreviewString=[PreviewString, TmpString];
        end
    end
    CorrSeedFile=get(handles.CorrSeedListbox, 'String');
    if ~isempty(CorrSeedFile)
        CorrSeedData=load(CorrSeedFile{1});
        [CorrPath, CorrName, CorrExt]=fileparts(CorrSeedFile{1});
        NumOfCorrSeed=size(CorrSeedData, 2);
        for j=1:NumOfCorrSeed
            ColumnNames=[ColumnNames; {sprintf('Corr Seed - %s(%d)', CorrName, j)}];
        end
        CorrSeedData=cellfun(@(n) num2str(n), num2cell(CorrSeedData),...
            'UniformOutput', false);
        TmpString=cell(L, NumOfCorrSeed);
        TmpString(1:size(CorrSeedData, 1), :)=CorrSeedData;
        PreviewString=[PreviewString, TmpString];        
    end
else
    PreviewString=cell(10, 10);
    ColumnNames=cell(10, 1);
end
set(handles.GLMTable, 'Data', PreviewString);
set(handles.GLMTable, 'ColumnName', ColumnNames);
Width=cell(1, numel(ColumnNames));
Width=cellfun(@(w) 150, Width, 'UniformOutput', false);
set(handles.GLMTable, 'ColumnWidth', Width);

function FileName=GetFileName(FullPath)
[Path, Name, Ext]=fileparts(FullPath);
FileName=[Name, Ext];

%% Get X and Y
function [Y_Volume, C_Volume, C_Text, R_Text, NiiHeader, NumOfSamples]=GetVolume(handles)
SampleCells=handles.SampleCells;
CovImageCells=handles.CovImageCells;
TextCells=get(handles.CovTextListbox, 'String');
RCells=get(handles.CorrSeedListbox, 'String');

NumOfSamples=cellfun(@numel, SampleCells);
NumSum=sum(NumOfSamples);

NiiHeader=[];
% Y Volume
Y_Cell=cell(1, 1, 1, NumSum);
for i=1:numel(NumOfSamples)
    % Obtain Begin and End Point
    if i==1
        B=1;
    else
        B=sum(NumOfSamples(1:i-1))+1;
    end
    E=sum(NumOfSamples(1:i));
    
    % Put Data into Cell
    OneGroup=SampleCells{i};
    if isempty(NiiHeader)
        NiiHeader=spm_vol(OneGroup{1});
    end
    fprintf('Reading Group-%d Image Sample\n', i);
    cellfun(@(s) fprintf('\t%s\n', s), OneGroup);
    Y_Cell(1, 1, 1, B:E)=cellfun(@(s) y_ReadRPI(s), OneGroup, 'UniformOutput', false);
end
Y_Volume=cell2mat(Y_Cell);

% Image Covariate
C_Volume=[];
if ~isempty(CovImageCells)
    if ~all(NumOfSamples==cellfun(@numel, CovImageCells))
        errordlg('Invalid Dimension of Covariate Image');
    end
    C_Cell=cell(1, 1, 1, NumSum);
    for i=1:numel(NumOfSamples)
        % Obtain Begin and End Point
        if i==1
            B=1;
        else
            B=NumOfSamples(i-1)+1;
        end
        E=sum(NumOfSamples(1:i));
    
        % Put Data into Cell
        OneImageCov=CovImageCells{i};
        fprintf('Reading Group-%d Image Covarites\n', i);
        cellfun(@(s) fprintf('\t%s\n', s), OneImageCov);
        C_Cell(1, 1, 1, B:E)=cellfun(@(s) y_ReadRPI(s), OneImageCov, 'UniformOutput', false);
    end
    C_Volume=cell2mat(C_Cell);  
end

% Text Covariate 
C_Text=[];
if ~isempty(TextCells)
    fprintf('Reading Text Covarites\n');
    cellfun(@(s) fprintf('\t%s\n', s), TextCells);
    C_Text=cellfun(@(s) load(s), TextCells, 'UniformOutput', false);
    if ~all(NumOfSamples==cellfun(@(d) size(d, 1), C_Text)')
        errordlg('Invalid Dimension of Covariate Text');
    end
    C_Text=cell2mat(C_Text);
end

R_Text=[];
if ~isempty(RCells)
    fprintf('Reading Correlation Seed\n');
    cellfun(@(s) fprintf('\t%s\n', s), RCells);
    R_Text=cellfun(@(s) load(s), RCells, 'UniformOutput', false);
    if ~all(NumOfSamples==cellfun(@(d) size(d, 1), R_Text))
        errordlg('Invalid Dimension of Correlated Seed');
    end
    R_Text=cell2mat(R_Text);
end

function [S, R]=Regression(Y_Volume, Regressor, C_Text, C_Volume, Contrast, TFR_Flag)
S=zeros(size(Y_Volume, 1), 1);
R=zeros(size(Y_Volume));

Xinit=Regressor;
if ~isempty(C_Text)
    Xinit=[Xinit, C_Text];
end

fprintf('Begin Regression');
for i=1:size(Y_Volume, 1)
    Y=Y_Volume(i, :)';
    
    if ~isempty(C_Volume)
        X=[Xinit, C_Volume(i, :)'];
    else
        X=Xinit;
    end
    [b, r, SSE, SSR, T, TF_ForContrast]=regress_ss(Y, X, Contrast, TFR_Flag);
    S(i, 1)=TF_ForContrast;
    R(i, :)=r;
    if ~rem(i, floor(size(Y_Volume, 1)/10));
        fprintf('.');
    end
end
fprintf('\n')
fprintf('Done\n');
S(isnan(S))=0;
R(isnan(R))=0;

function [b,r,SSE,SSR, T, TFR_ForContrast] = regress_ss(y,X,Contrast,TFR_Flag)
% [b,r,SSE,SSR, T, TF_ForContrast] = y_regress_ss(y,X,Contrast,TF_Flag)
% Perform regression.
% Revised from MATLAB's regress in order to speed up the calculation.
% Input:
%   y - Independent variable.
%   X - Dependent variable.
%   Contrast [optional] - Contrast for T-test for F-test. 1*ncolX matrix.
%   TF_Flag [optional] - 'T' or 'F'. Specify if T-test or F-test need to be performed for the contrast
% Output:
%   b - beta of regression model.
%   r - residual.
%   SSE - The sum of squares of error.
%   SSR - The sum of squares of regression.
%   T - T value for each beta.
%   TF_ForContrast - T or F value (depends on TF_Flag) for the contrast.
%   
%___________________________________________________________________________
% Written by YAN Chao-Gan 100317.
% The Nathan Kline Institute for Psychiatric Research, 140 Old Orangeburg Road, Orangeburg, NY 10962; 
% Child Mind Institute, 445 Park Avenue, New York, NY 10022; 
% The Phyllis Green and Randolph Cowen Institute for Pediatric Neuroscience, New York University Child Study Center, New York, NY 10016
% State Key Laboratory of Cognitive Neuroscience and Learning, Beijing Normal University, China, 100875
% ycg.yan@gmail.com
% Revised by YAN Chao-Gan 120519. Also output T value for each beta. Referenced from regstats.m
% Revised by YAN Chao-Gan 121220. Also support T-test or F-test for a given contrast.
% Revised y_regress_ss in DPABI by Sandy

[n,ncolX] = size(X);
[Q,R,perm] = qr(X,0);
p = sum(abs(diag(R)) > max(n,ncolX)*eps(R(1)));
if p < ncolX,
    R = R(1:p,1:p);
    Q = Q(:,1:p);
    perm = perm(1:p);
end
b = zeros(ncolX,1);
b(perm) = R \ (Q'*y);

if nargout >= 2
    yhat = X*b;                     % Predicted responses at each data point.
    r = y-yhat;                     % Residuals.
    if nargout >= 3
        SSE=sum(r.^2);
        if nargout >= 4
            SSR=sum((yhat-mean(y)).^2);
            
            if nargout >= 5
                %Also output T value for each beta. Referenced from regstats.m
                [Q,R] = qr(X,0);
                ri = R\eye(ncolX);
                T = b./sqrt(diag(ri*ri' * (SSE/(n-ncolX))));
                
                if nargout >= 6
                    %YAN Chao-Gan 121220. Also support T-test or F-test for a given contrast.
                    % Have contrast
                    if strcmpi(TFR_Flag,'T') || strcmpi(TFR_Flag,'R') %Add R Flag by Sandy
                        std_e = sqrt(SSE/(n-ncolX));        % Standard deviation of the noise
                        d = sqrt(Contrast*(X'*X)^(-1)*Contrast');
                        TFR_ForContrast = (Contrast*b)./(std_e*d);           % T-test
                        if strcmpi(TFR_Flag, 'R')
                            Df_E = size(X, 1) - size(Contrast,2);
                            TFR_ForContrast=TFR_ForContrast./(sqrt(Df_E+TFR_ForContrast.*TFR_ForContrast));
                        end
                    elseif strcmpi(TFR_Flag,'F')
                        X0 = X(:,~Contrast);
                        ncolX0 = size(X0,2);
                        if ncolX0>0
                            b0 = (X0'*X0)^(-1)*X0'*y; % Regression coefficients (restricted model)
                            r0 = y-X0*b0;
                            SSE0 = sum(r0.^2); % Estimate of the residual sum-of-square of the restricted model (SSE0)
                        else
                            SSE0 = sum(y.^2,1);
                        end
                        TFR_ForContrast = ((SSE0-SSE)/(ncolX-ncolX0))./(SSE/(n-ncolX)); % F-Test
                        
                    end
                end
                
            end
            
        end
    end
end

function [dLh, FWHM]=Smoothest(R_Volume, MaskData, DOF, Vox)
%Revised by Sandy from DPABI y_Smoothest.m
fprintf('Begin Smoothest\n')
R_Volume=single(R_Volume);
[N1, N2, N3, N4]=size(R_Volume);

if N4>2
    R_Volume=(R_Volume-repmat(mean(R_Volume,4),[1,1,1, N4]))./repmat(std(R_Volume,0,4),[1,1,1, N4]);%Zero mean and one std
    R_Volume(isnan(R_Volume))=0;
end

SSminus=[0 0 0];
S2=[0 0 0];

N=0;
for x=2:N1
    for y=2:N2
        for z=2:N3
            if MaskData(x, y, z) && MaskData(x-1, y, z) && MaskData(x, y-1, z) && MaskData(x, y, z-1)
                N=N+1;
                for t=1:N4
                    SSminus(1) = SSminus(1) + R_Volume(x, y, z, t) * R_Volume(x-1, y, z, t);
                    SSminus(2) = SSminus(2) + R_Volume(x, y, z, t) * R_Volume(x, y-1, z, t);
                    SSminus(3) = SSminus(3) + R_Volume(x, y, z, t) * R_Volume(x, y, z-1, t);
                    
                    S2(1) = S2(1) + 0.5 * ((R_Volume(x, y, z, t)^2) + (R_Volume(x-1, y, z, t)^2));
                    S2(2) = S2(2) + 0.5 * ((R_Volume(x, y, z, t)^2) + (R_Volume(x, y-1, z, t)^2));
                    S2(3) = S2(3) + 0.5 * ((R_Volume(x, y, z, t)^2) + (R_Volume(x, y, z-1, t)^2));
                end
            end
        end
    end
    fprintf('.');
end

if SSminus(1)>0.99999999*S2(1)
    SSminus(1)=0.99999999*S2(1);
    warning('possibly biased smootheness in X');
end
if SSminus(2)>0.99999999*S2(2)
    SSminus(2)=0.99999999*S2(2);
    warning('possibly biased smootheness in Y');
end
if SSminus(3)>0.99999999*S2(3)
    SSminus(3)=0.99999999*S2(3);
    warning('possibly biased smootheness in Z');
end

sigmasq(1) = -1 / (4 * log(abs(SSminus(1)/S2(1))));
sigmasq(2) = -1 / (4 * log(abs(SSminus(2)/S2(2))));
sigmasq(3) = -1 / (4 * log(abs(SSminus(3)/S2(3))));

dLh=((sigmasq(1)*sigmasq(2)*sigmasq(3))^-0.5)*(8^-0.5);

if N4 > 1
    fprintf('DLH %f voxels^-3 before correcting for temporal DOF\n',dLh);
    
    lut(6)   = 1.5423138; lut(7)   = 1.3757105; lut(8)   = 1.2842680;
    lut(9)   = 1.2272151; lut(10)  = 1.1885232; lut(11)  = 1.1606988;
    lut(12)  = 1.1398000; lut(13)  = 1.1235677; lut(14)  = 1.1106196;
    lut(15)  = 1.1000651; lut(16)  = 1.0913060; lut(17)  = 1.0839261;
    lut(18)  = 1.0776276; lut(19)  = 1.0721920; lut(20)  = 1.0674553;
    lut(21)  = 1.0632924; lut(26)  = 1.0483053; lut(31)  = 1.0390117;
    lut(41)  = 1.0281339; lut(51)  = 1.0219834; lut(61)  = 1.0180339;
    lut(71)  = 1.0152850; lut(81)  = 1.0132621; lut(91)  = 1.0117115;
    lut(101) = 1.0104851; lut(151) = 1.0068808; lut(201) = 1.0051200;
    lut(301) = 1.0033865; lut(501) = 1.0020191;
    
    y = lut(lut~=0);
    x = find(lut~=0);
    xi=[1:501];
    lut_interpolated=interp1(x,y,xi,'linear');
    
    if (DOF < 6)
        dLh=dLh * 1.1;
    elseif (DOF>500)
        dLh=dLh * (1.0321/DOF +1)^0.5;
    else
        retval=(lut_interpolated(floor(DOF)+1)-lut_interpolated(floor(DOF)))*(floor(DOF)+1-floor(DOF)) + ...
            lut_interpolated(floor(DOF)+1);
        dLh=dLh * retval^0.5;
    end
    
end

FWHM(1) =  sqrt(8 * log(2) * sigmasq(1));
FWHM(2) =  sqrt(8 * log(2) * sigmasq(2));
FWHM(3) =  sqrt(8 * log(2) * sigmasq(3));

resels = FWHM(1)*FWHM(2)*FWHM(3);
fprintf('\nFWHMx = %f voxels\nFWHMy = %f voxels\nFWHMz = %f voxels\n',FWHM(1),FWHM(2),FWHM(3));
FWHM=FWHM.*Vox;
fprintf('FWHMx = %f mm\nFWHMy = %f mm\nFWHMz = %f mm\n',FWHM(1),FWHM(2),FWHM(3));
nVoxels=length(find(MaskData));
fprintf('DLH = %f\nVOLUME = %d\nRESELS = %f\n',dLh,nVoxels,resels);
