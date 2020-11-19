function varargout = SeeCAT_FuncConnect(varargin)
% SEECAT_FUNCCONNECT MATLAB code for SeeCAT_FuncConnect.fig
%      SEECAT_FUNCCONNECT, by itself, creates a new SEECAT_FUNCCONNECT or raises the existing
%      singleton*.
%
%      H = SEECAT_FUNCCONNECT returns the handle to a new SEECAT_FUNCCONNECT or the handle to
%      the existing singleton*.
%
%      SEECAT_FUNCCONNECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEECAT_FUNCCONNECT.M with the given input arguments.
%
%      SEECAT_FUNCCONNECT('Property','Value',...) creates a new SEECAT_FUNCCONNECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SeeCAT_FuncConnect_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SeeCAT_FuncConnect_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SeeCAT_FuncConnect

% Last Modified by GUIDE v2.5 05-Jul-2017 14:49:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SeeCAT_FuncConnect_OpeningFcn, ...
    'gui_OutputFcn',  @SeeCAT_FuncConnect_OutputFcn, ...
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


% --- Executes just before SeeCAT_FuncConnect is made visible.
function SeeCAT_FuncConnect_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SeeCAT_FuncConnect (see VARARGIN)

% Choose default command line output for SeeCAT_FuncConnect
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

movegui(hObject,'center');
% UIWAIT makes SeeCAT_FuncConnect wait for user response (see UIRESUME)
% uiwait(handles.SeeCAT_FC_figure);
set(handles.WorkDir_edit,'String',pwd);
ProgramPath = fileparts(which('SeeCAT.m'));
if ispc
    filesep = '\';
else
    filesep = '/';
end
set(handles.Mask_edit,'String',[ProgramPath,filesep,'templates',filesep,'GreyMask_02_61x73x61.nii']);



% --- Outputs from this function are returned to the command line.
function varargout = SeeCAT_FuncConnect_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function WorkDir_edit_Callback(hObject, eventdata, handles)
% hObject    handle to WorkDir_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WorkDir_edit as text
%        str2double(get(hObject,'String')) returns contents of WorkDir_edit as a double


% --- Executes during object creation, after setting all properties.
function WorkDir_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WorkDir_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in WorkDir_pushbutton.
function WorkDir_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to WorkDir_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder_name = uigetdir(pwd,'Select Working Directory');
if isequal(folder_name,0)
    return;
else
    set(handles.WorkDir_edit,'String',folder_name);
end



function StartFolder_edit_Callback(hObject, eventdata, handles)
% hObject    handle to StartFolder_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StartFolder_edit as text
%        str2double(get(hObject,'String')) returns contents of StartFolder_edit as a double


% --- Executes during object creation, after setting all properties.
function StartFolder_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StartFolder_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Seed_listbox.
function Seed_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to Seed_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Seed_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Seed_listbox
persistent chk
if isempty(chk)
    chk = 1;
    pause(0.5);
    if chk == 1
        chk = [];
    end
else
    chk = [];
    list_string = get(hObject,'String');
    list_pos = get(hObject,'Value');
    list_string(list_pos) = [];
    if isempty(list_string)
        list_string{1} = [];
    end
    if list_pos > length(list_string)
        set(hObject,'Value',list_pos - 1);
    end
    set(hObject,'String',list_string);
end

% --- Executes during object creation, after setting all properties.
function Seed_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Seed_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AddSpheres_pushbutton.
function AddSpheres_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to AddSpheres_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SphereMat = AddSpheres;
SphereList = num2cell(SphereMat,2);
SeedList = get(handles.Seed_listbox,'String');
nSeed = length(SeedList);
if nSeed == 1 && isempty(SeedList{1})
    nSeed = 0;
end
for i = nSeed+1:nSeed+length(SphereList)
    SeedList{i,1} = ['Sphere: ',num2str(SphereList{i-nSeed})];
end
set(handles.Seed_listbox,'String',SeedList);



% --- Executes on button press in AddClusters_pushbutton.
function AddClusters_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to AddClusters_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Name, Path] = uigetfile({'*.img;*.nii;*.nii.gz','Brain Image Files (*.img;*.nii;*.nii.gz)';'*.*', 'All Files (*.*)';},...
    'Pick the Masks for ROI', 'MultiSelect','on');
if isnumeric(Name)
    return
end
SeedList = get(handles.Seed_listbox,'String');
nSeed = length(SeedList);
if nSeed == 1 && isempty(SeedList{1})
    nSeed = 0;
end
if ischar(Name)
    Name={Name};
end
for i = nSeed+1:nSeed+length(Name)
    SeedList{i,1} = ['Mask: ',Path,Name{i-nSeed}];
end
set(handles.Seed_listbox,'String',SeedList);




function Mask_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Mask_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Mask_edit as text
%        str2double(get(hObject,'String')) returns contents of Mask_edit as a double


% --- Executes during object creation, after setting all properties.
function Mask_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mask_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Mask_pushbutton.
function Mask_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Mask_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Name, Path] = uigetfile({'*.img;*.nii;*.nii.gz','Brain Image Files (*.img;*.nii;*.nii.gz)';'*.*', 'All Files (*.*)';},...
    'Pick the Masks for Calculation');
if isnumeric(Name)
    return
end
set(handles.Mask_edit,'String',[Path,Name]);


function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Run_pushbutton.
function Run_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Run_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'Enable','off');
set(hObject,'String','running...');
tic
if get(handles.Parallel_checkbox,'Value') == 1
    if exist('gcp.m','file')
        try
            gcp;
        end
    elseif  matlabpool('size') == 0
        try
            matlabpool;
        end
    end
end
if ispc
    filesep = '\';
else
    filesep = '/';
end
WorkDir = get(handles.WorkDir_edit,'String');
StartFolder = get(handles.StartFolder_edit,'String');
SeedList = get(handles.Seed_listbox,'String');
nSeed = length(SeedList);
MaskName = get(handles.Mask_edit,'String');

Log_Text = datestr(now);
Log_Text = sprintf('%s\n%s\n',Log_Text,['Working directory: ',WorkDir]);
Log_Text = sprintf('%s%s\n\n',Log_Text,['Start folder: ',StartFolder]);
Log_Text = sprintf('%s%s\n',Log_Text,'Seed list: ');
for i = 1:nSeed
    Log_Text = sprintf('%s%s\n',Log_Text,['  ',SeedList{i}]);
end
Log_Text = sprintf('%s%s\n',Log_Text,'');
Log_Text = sprintf('%s%s\n',Log_Text,['Calculation mask: ',MaskName]);
fid = fopen([WorkDir,filesep,'Log_FC_',datestr(now,'yyyymmddHHMMSS'),'.txt'],'at+');
fprintf(fid,'%s',Log_Text);
fclose(fid);


fprintf('Calculating Functional Connectivity Started.\n');

hdr_mask = spm_vol(MaskName);
[vol_mask,XYZ] = spm_read_vols(hdr_mask);

vol_seed = cell(nSeed,1);
for i = 1:nSeed
    switch SeedList{i}(1)
        case 'S'
            CoordMat = str2num(SeedList{i}(8:end));
            distance = pdist2(CoordMat(1:3),XYZ');
            vol_seed{i} = zeros(hdr_mask.dim);
            vol_seed{i}(distance<=CoordMat(4)) = 1;
        case 'M'
            hdr_tmp = spm_vol(SeedList{i}(7:end));
            vol_seed{i} = spm_read_vols(hdr_tmp);
    end
end

Sublist = dir([WorkDir,filesep,StartFolder]);
if strcmpi(Sublist(3).name,'.DS_Store')
    Sublist(1:3) = [];
else
    Sublist(1:2) = [];
end
mask_ind = reshape(vol_mask>0,1,[]);
mkdir([WorkDir,filesep,'Results',filesep,'FC_',StartFolder]);
if get(handles.Parallel_checkbox,'Value') == 1
    parfor i = 1:length(Sublist)
        cd([WorkDir,filesep,StartFolder,filesep,Sublist(i).name]);
        Filename = dir('*.nii');
        Nii = nifti(Filename.name);
        volCourse = reshape(double(Nii.dat),[Nii.dat.dim(1,1) * Nii.dat.dim(1,2) * Nii.dat.dim(1,3),Nii.dat.dim(1,4)])';
        maskCourse = volCourse(:,mask_ind);
        seedCourse = zeros(Nii.dat.dim(1,4),nSeed);
        for j = 1:nSeed
            seedCourse(:,j) = mean(volCourse(:,reshape(vol_seed{j}>0,1,[])),2);
        end
        seedCourse = seedCourse - repmat(mean(seedCourse),Nii.dat.dim(1,4),1);
        seedCourse = seedCourse./repmat(std(seedCourse,0,1),Nii.dat.dim(1,4),1);
        maskCourse = maskCourse - repmat(mean(maskCourse),Nii.dat.dim(1,4),1);
        maskCourse = maskCourse./repmat(std(maskCourse,0,1),Nii.dat.dim(1,4),1);
        r = seedCourse' * maskCourse ./(Nii.dat.dim(1,4)-1);
        z = FisherTrans(r);
        cd([WorkDir,filesep,'Results',filesep,'FC_',StartFolder]);
        for j = 1:nSeed
            hdr_fc = hdr_mask;
            Num = ['00',num2str(j)];
            Num = Num(end-1:end);
            hdr_fc.fname = ['FC_Seed',Num,'_',Sublist(i).name,'.nii'];
            hdr_fc.dt(1) = 16;
            vol_fc = zeros(hdr_fc.dim);
            vol_fc(mask_ind) = r(j,:);
            spm_write_vol(hdr_fc,vol_fc);
            hdr_zfc = hdr_fc;
            hdr_zfc.fname = ['z',hdr_zfc.fname];
            vol_zfc = zeros(hdr_zfc.dim);
            vol_zfc(mask_ind) = z(j,:);
            spm_write_vol(hdr_zfc,vol_zfc);
        end
    end
else
    for i = 1:length(Sublist)
        cd([WorkDir,filesep,StartFolder,filesep,Sublist(i).name]);
        Filename = dir('*.nii');
        Nii = nifti(Filename.name);
        volCourse = reshape(double(Nii.dat),[Nii.dat.dim(1,1) * Nii.dat.dim(1,2) * Nii.dat.dim(1,3),Nii.dat.dim(1,4)])';
        maskCourse = volCourse(:,mask_ind);
        seedCourse = zeros(Nii.dat.dim(1,4),nSeed);
        for j = 1:nSeed
            seedCourse(:,j) = mean(volCourse(:,reshape(vol_seed{j}>0,1,[])),2);
        end
        seedCourse = seedCourse - repmat(mean(seedCourse),Nii.dat.dim(1,4),1);
        seedCourse = seedCourse./repmat(std(seedCourse,0,1),Nii.dat.dim(1,4),1);
        maskCourse = maskCourse - repmat(mean(maskCourse),Nii.dat.dim(1,4),1);
        maskCourse = maskCourse./repmat(std(maskCourse,0,1),Nii.dat.dim(1,4),1);
        r = seedCourse' * maskCourse ./(Nii.dat.dim(1,4)-1);
        z = FisherTrans(r);
        cd([WorkDir,filesep,'Results',filesep,'FC_',StartFolder]);
        for j = 1:nSeed
            hdr_fc = hdr_mask;
            Num = ['00',num2str(j)];
            Num = Num(end-1:end);
            hdr_fc.fname = ['FC_Seed',Num,'_',Sublist(i).name,'.nii'];
            hdr_fc.dt(1) = 16;
            vol_fc = zeros(hdr_fc.dim);
            vol_fc(mask_ind) = r(j,:);
            spm_write_vol(hdr_fc,vol_fc);
            hdr_zfc = hdr_fc;
            hdr_zfc.fname = ['z',hdr_zfc.fname];
            vol_zfc = zeros(hdr_zfc.dim);
            vol_zfc(mask_ind) = z(j,:);
            spm_write_vol(hdr_zfc,vol_zfc);
        end
    end
end
cd([WorkDir,filesep,'Results',filesep,'FC_',StartFolder]);
for i = 1:nSeed
    Num = ['00',num2str(i)];
    Num = Num(end-1:end);
    movefile(['zFC_Seed',Num,'*.nii'],[WorkDir,filesep,'Results',filesep,'FC_',StartFolder,filesep,'zFC_Seed',Num]);
    movefile(['FC_Seed',Num,'*.nii'],[WorkDir,filesep,'Results',filesep,'FC_',StartFolder,filesep,'FC_Seed',Num]);
end
cd(WorkDir);
fprintf('Calculating Functional Connectivity Finished.\n');
toc
set(hObject,'Enable','on');
set(hObject,'String','RUN');


function z = FisherTrans(r)
z = 0.5*log((1+r)./(1-r));


% --- Executes on button press in Quit_pushbutton.
function Quit_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Quit_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.SeeCAT_FC_figure);


% --- Executes on button press in Parallel_checkbox.
function Parallel_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Parallel_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Parallel_checkbox
