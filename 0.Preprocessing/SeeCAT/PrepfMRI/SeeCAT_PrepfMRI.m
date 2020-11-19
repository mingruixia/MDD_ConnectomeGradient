function varargout = SeeCAT_PrepfMRI(varargin)
% SEECAT_PREPFMRI MATLAB code for SeeCAT_PrepfMRI.fig
%      SEECAT_PREPFMRI, by itself, creates a new SEECAT_PREPFMRI or raises the existing
%      singleton*.
%
%      H = SEECAT_PREPFMRI returns the handle to a new SEECAT_PREPFMRI or the handle to
%      the existing singleton*.
%
%      SEECAT_PREPFMRI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEECAT_PREPFMRI.M with the given input arguments.
%
%      SEECAT_PREPFMRI('Property','Value',...) creates a new SEECAT_PREPFMRI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SeeCAT_PrepfMRI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SeeCAT_PrepfMRI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SeeCAT_PrepfMRI

% Last Modified by GUIDE v2.5 04-Jul-2017 10:56:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SeeCAT_PrepfMRI_OpeningFcn, ...
    'gui_OutputFcn',  @SeeCAT_PrepfMRI_OutputFcn, ...
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


% --- Executes just before SeeCAT_PrepfMRI is made visible.
function SeeCAT_PrepfMRI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SeeCAT_PrepfMRI (see VARARGIN)

% Choose default command line output for SeeCAT_PrepfMRI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

movegui(hObject,'center');

% UIWAIT makes SeeCAT_PrepfMRI wait for user response (see UIRESUME)
% uiwait(handles.SeeCAT_PrefMRI_figure);
set(handles.WorkDir_edit,'String',pwd);


if ispc
    filesep = '\';
else
    filesep = '/';
end
list = dir([get(handles.WorkDir_edit,'String'),filesep,get(handles.StartFolder_edit,'String')]);
if ~isempty(list)
    if strcmpi(list(3).name,'.DS_Store')
        list(1:3) = [];
    else
        list(1:2) = [];
    end
    list_string = cell(1,length(list));
    for i = 1:length(list);
        list_string{i} = list(i).name;
    end
    set(handles.SubList_listbox,'String',list_string);
end


% --- Outputs from this function are returned to the command line.
function varargout = SeeCAT_PrepfMRI_OutputFcn(hObject, eventdata, handles)
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
list = dir([folder_name,filesep,get(handles.StartFolder_edit,'String')]);
if ~isempty(list)
    if strcmpi(list(3).name,'.DS_Store')
        list(1:3) = [];
    else
        list(1:2) = [];
    end
    list_string = cell(1,length(list));
    for i = 1:length(list);
        list_string{i} = list(i).name;
    end
    set(handles.SubList_listbox,'String',list_string);
end


function StartFolder_edit_Callback(hObject, eventdata, handles)
% hObject    handle to StartFolder_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StartFolder_edit as text
%        str2double(get(hObject,'String')) returns contents of StartFolder_edit as a double

if ispc
    filesep = '\';
else
    filesep = '/';
end
list = dir([get(handles.WorkDir_edit,'String'),filesep,get(hObject,'String')]);
if isempty(list)
    warndlg(['There is no ',get(hObject,'String'),' folder under working directory'],'Warning');
    list_string = [];
    set(handles.SubList_listbox,'String',list_string);
else
    if strcmpi(list(3).name,'.DS_Store')
        list(1:3) = [];
    else
        list(1:2) = [];
    end
    list_string = cell(1,length(list));
    for i = 1:length(list);
        list_string{i} = list(i).name;
    end
    set(handles.SubList_listbox,'String',list_string);
end


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


% --- Executes on selection change in SubList_listbox.
function SubList_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to SubList_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SubList_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SubList_listbox
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
function SubList_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SubList_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TR_edit_Callback(hObject, eventdata, handles)
% hObject    handle to TR_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TR_edit as text
%        str2double(get(hObject,'String')) returns contents of TR_edit as a double


% --- Executes during object creation, after setting all properties.
function TR_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TR_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Run_pushbutton.
function Run_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Run_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% set(hObject,'BackgroundColor',[153,255,255]/255);
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
Sublist = get(handles.SubList_listbox,'String');
SubNum = length(Sublist);
ProgramPath = fileparts(which('SeeCAT_PrepfMRI.m'));
DataSegNum = 10;

Log_Text = datestr(now);
Log_Text = sprintf('%s\n%s\n',Log_Text,['Working directory: ',WorkDir]);
Log_Text = sprintf('%s%s\n',Log_Text,['Start folder: ',StartFolder]);
Log_Text = sprintf('%s%s\n',Log_Text,['TR = ',get(handles.TR_edit,'String')]);
if get(handles.DICOMtoNIFTI_checkbox,'Value') == 1
    Log_Text = sprintf('%s%s\n',Log_Text,'DICOM to NIFTI');
end
if get(handles.RemoveTPS_checkbox,'Value') == 1
    Log_Text = sprintf('%s%s\n',Log_Text,['Remove first ', get(handles.RemoveTPS_edit,'String'),' time points']);
end
if get(handles.SliceTiming_checkbox,'Value') == 1
    slicetiming_text = get(handles.SliceTiming_popupmenu,'String');
    Log_Text = sprintf('%s%s\n',Log_Text,['Slice timing: ', slicetiming_text{get(handles.SliceTiming_popupmenu,'Value')}]);
end
if get(handles.Realign_checkbox,'Value') == 1
    Log_Text = sprintf('%s%s\n',Log_Text,'Realign');
end
if get(handles.Normalize_checkbox,'Value') == 1
    normalize_text = get(handles.Normalize_popupmenu,'String');
    Log_Text = sprintf('%s%s\n',Log_Text,['Normalize: using ',normalize_text{get(handles.Normalize_popupmenu,'Value')},'; Voxel size = [',get(handles.Normalize_edit,'String'),']']);
end
if get(handles.Smooth_checkbox,'Value') == 1
    Log_Text = sprintf('%s%s\n',Log_Text,['Smooth: FWHM = ',get(handles.Smooth_edit,'String'),']']);
end
if get(handles.Detrend_checkbox,'Value') == 1
    Log_Text = sprintf('%s%s\n',Log_Text,'Detrend');
end
if get(handles.Regress_checkbox,'Value') == 1
    Log_Text = sprintf('%s%s',Log_Text,'Regress Covariates: ');
    if get(handles.RegHM_checkbox,'Value') == 1
        reghm_text = get(handles.RegHM_popupmenu,'String');
        Log_Text = sprintf('%s%s',Log_Text,[reghm_text{get(handles.RegHM_popupmenu,'Value')},'; ']);
    end
    if get(handles.RegWM_checkbox,'Value') == 1
        Log_Text = sprintf('%s%s',Log_Text,'White Matter; ');
    end
    if get(handles.RegCSF_checkbox,'Value') == 1
        Log_Text = sprintf('%s%s',Log_Text,'Cerebrospinal Fluid; ');
    end
    if get(handles.RegGS_checkbox,'Value') == 1
        Log_Text = sprintf('%s%s',Log_Text,'Global Signal. ');
    end
    Log_Text = sprintf('%s%s\n',Log_Text,' ');
end
if get(handles.Filter_checkbox,'Value') == 1
    Log_Text = sprintf('%s%s\n',Log_Text,['Filter: ',get(handles.FilterLow_edit,'String'),'~',get(handles.FilterHigh_edit,'String'),' Hz']);
end

if get(handles.Scrubbing_checkbox,'Value') == 1
    scrubbing_text = get(handles.Scrubbing_popupmenu,'String');
    Log_Text = sprintf('%s%s\n',Log_Text,['Scrubbing: Threshold ',...
        get(handles.ScrubbingThr_edit,'String'),' Method ',...
        scrubbing_text{get(handles.Scrubbing_popupmenu,'Value')},...
        ' Before tps ',get(handles.ScrubbingBefore_edit,'String'),...
        ' After tps ',get(handles.ScrubbingAfter_edit,'String')]);
end

fid = fopen([WorkDir,filesep,'Log_Preprocess_',datestr(now,'yyyymmddHHMMSS'),'.txt'],'at+');
fprintf(fid,'%s',Log_Text);
fclose(fid);


% DICOM to NIFTI
if get(handles.DICOMtoNIFTI_checkbox,'Value') == 1
    % Convert Functional Images
    fprintf('Converting Functional Images Started.\n');
    if get(handles.Parallel_checkbox,'Value') == 1
        parfor i = 1:SubNum
%             delete([WorkDir,filesep,'FunRaw',filesep,'*.nii']);
%             delete([WorkDir,filesep,'FunRaw',filesep,'*.gz']);
%             delete([WorkDir,filesep,'FunRaw',filesep,'*.img']);
%             delete([WorkDir,filesep,'FunRaw',filesep,'*.hdr']);
            DCMfile = dir([WorkDir,filesep,'FunRaw',filesep,Sublist{i}]);
            if strcmpi(DCMfile(3).name,'.DS_Store')
                StartInd=4;
            else
                StartInd=3;
            end
            Filename = [WorkDir,filesep,'FunRaw',filesep,Sublist{i},filesep,DCMfile(StartInd).name];
            OutputDir = [WorkDir,filesep,'FunImg',filesep,Sublist{i}];
            mkdir(OutputDir);
            fprintf(['Converting Functional Images:',Sublist{i},' Started.\n']);
            dcm2nii(Filename, OutputDir);
            fprintf(['Converting Functional Images:',Sublist{i},' OK.\n']);
        end
    else
        for i = 1:SubNum
%             delete([WorkDir,filesep,'FunRaw',filesep,'*.nii']);
%             delete([WorkDir,filesep,'FunRaw',filesep,'*.gz']);
%             delete([WorkDir,filesep,'FunRaw',filesep,'*.img']);
%             delete([WorkDir,filesep,'FunRaw',filesep,'*.hdr']);
            DCMfile = dir([WorkDir,filesep,'FunRaw',filesep,Sublist{i}]);
            if strcmpi(DCMfile(3).name,'.DS_Store')
                StartInd=4;
            else
                StartInd=3;
            end
            Filename = [WorkDir,filesep,'FunRaw',filesep,Sublist{i},filesep,DCMfile(StartInd).name];
            OutputDir = [WorkDir,filesep,'FunImg',filesep,Sublist{i}];
            mkdir(OutputDir);
            fprintf(['Converting Functional Images:',Sublist{i},' Started.\n']);
            dcm2nii(Filename, OutputDir);
            fprintf(['Converting Functional Images:',Sublist{i},' OK.\n']);
        end
    end
    fprintf('Converting Functional Images Finished.\n');
    StartFolder = 'FunImg';
    
    % Convert T1 Images
    if ~isempty(dir([WorkDir,filesep,'T1Raw'])) && isempty(dir([WorkDir,filesep,'T1Img']))
        fprintf('Converting T1 Images Started.\n');
        if get(handles.Parallel_checkbox,'Value') == 1
            parfor i = 1:SubNum
                DCMfile = dir([WorkDir,filesep,'T1Raw',filesep,Sublist{i}]);
                if strcmpi(DCMfile(3).name,'.DS_Store')
                    StartInd=4;
                else
                    StartInd=3;
                end
                Filename = [WorkDir,filesep,'T1Raw',filesep,Sublist{i},filesep,DCMfile(StartInd).name];
                OutputDir = [WorkDir,filesep,'T1Img',filesep,Sublist{i}];
                mkdir(OutputDir);
                fprintf(['Converting T1 Images:',Sublist{i},' Started.\n']);
                dcm2nii(Filename, OutputDir);
                fprintf(['Converting T1 Images:',Sublist{i},' OK.\n']);
            end
        else
            for i = 1:SubNum
                DCMfile = dir([WorkDir,filesep,'T1Raw',filesep,Sublist{i}]);
                if strcmpi(DCMfile(3).name,'.DS_Store')
                    StartInd=4;
                else
                    StartInd=3;
                end
                Filename = [WorkDir,filesep,'T1Raw',filesep,Sublist{i},filesep,DCMfile(StartInd).name];
                OutputDir = [WorkDir,filesep,'T1Img',filesep,Sublist{i}];
                mkdir(OutputDir);
                fprintf(['Converting T1 Images:',Sublist{i},' Started.\n']);
                dcm2nii(Filename, OutputDir);
                fprintf(['Converting T1 Images:',Sublist{i},' OK.\n']);
            end
        end
    end
    
    % Convert FieldMap Images
    if ~isempty(dir([WorkDir,filesep,'FieldRaw'])) && isempty(dir([WorkDir,filesep,'FieldImg']))
        fprintf('Converting FieldMap Images Started.\n');
        if get(handles.Parallel_checkbox,'Value') == 1
            parfor i = 1:SubNum
                DCMfile = dir([WorkDir,filesep,'FieldRaw',filesep,Sublist{i}]);
                if strcmpi(DCMfile(3).name,'.DS_Store')
                    StartInd=4;
                else
                    StartInd=3;
                end
                Filename = [WorkDir,filesep,'FieldRaw',filesep,Sublist{i},filesep,DCMfile(StartInd).name];
                OutputDir = [WorkDir,filesep,'FieldImg',filesep,Sublist{i}];
                mkdir(OutputDir);
                fprintf(['Converting FieldMap Images:',Sublist{i},' Started.\n']);
                dcm2nii(Filename, OutputDir);
                fprintf(['Converting FieldMap Images:',Sublist{i},' OK.\n']);
            end
        else
            for i = 1:SubNum
                DCMfile = dir([WorkDir,filesep,'FieldRaw',filesep,Sublist{i}]);
                if strcmpi(DCMfile(3).name,'.DS_Store')
                    StartInd=4;
                else
                    StartInd=3;
                end
                Filename = [WorkDir,filesep,'FieldRaw',filesep,Sublist{i},filesep,DCMfile(StartInd).name];
                OutputDir = [WorkDir,filesep,'FieldImg',filesep,Sublist{i}];
                mkdir(OutputDir);
                fprintf(['Converting FieldMap Images:',Sublist{i},' Started.\n']);
                dcm2nii(Filename, OutputDir);
                fprintf(['Converting FieldMap Images:',Sublist{i},' OK.\n']);
            end
        end
    end
    
    fprintf('Converting T1 Images Finished.\n');
    cd(WorkDir);
end

% Delete first n time points
if get(handles.RemoveTPS_checkbox,'Value') == 1
    DeleteNum = str2double(get(handles.RemoveTPS_edit,'String'));
    fprintf(['Removing First ',get(handles.RemoveTPS_edit,'String'),' Images Started.\n']);
    if get(handles.Parallel_checkbox,'Value') == 1
        parfor i=1:SubNum
            fprintf(['Removing First ',num2str(DeleteNum),' Time Points: ',Sublist{i},' Started.\n']);
            cd([WorkDir,filesep,StartFolder,filesep,Sublist{i}]);
            Filename = dir('*.nii');
            if ~isempty(Filename)
                if length(Filename) == 1
                    Nii  = nifti(Filename.name);
                    Data = Nii.dat(:,:,:,DeleteNum+1:end);
                    NewNii = Nii;
                    NewNii.dat.dim(1,4) = Nii.dat.dim(1,4) - DeleteNum;
                    create(NewNii);
                    NewNii.dat(:,:,:,:)=Data;
                else
                    for j = 1:DeleteNum
                        delete(Filename(j).name);
                    end
                end
            else
                Filename = dir('*.img');
                for j = 1:DeleteNum
                    delete(Filename(j).name);
                    delete([Filename(j).name(1:end-4),'.hdr']);
                end
            end
            fprintf(['Removing First ',num2str(DeleteNum),' Time Points: ',Sublist{i},' OK.\n']);
        end
    else
        for i=1:SubNum
            fprintf(['Removing First ',num2str(DeleteNum),' Time Points: ',Sublist{i},' Started.\n']);
            cd([WorkDir,filesep,StartFolder,filesep,Sublist{i}]);
            Filename = dir('*.nii');
            if ~isempty(Filename)
                if length(Filename) == 1
                    Nii  = nifti(Filename.name);
                    Data = Nii.dat(:,:,:,DeleteNum+1:end);
                    NewNii = Nii;
                    NewNii.dat.dim(1,4) = Nii.dat.dim(1,4) - DeleteNum;
                    create(NewNii);
                    NewNii.dat(:,:,:,:)=Data;
                else
                    for j = 1:DeleteNum
                        delete(Filename(j).name);
                    end
                end
            else
                Filename = dir('*.img');
                for j = 1:DeleteNum
                    delete(Filename(j).name);
                    delete([Filename(j).name(1:end-4),'.hdr']);
                end
            end
            fprintf(['Removing First ',num2str(DeleteNum),' Time Points: ',Sublist{i},' OK.\n']);
        end
    end
    fprintf(['Removing First ',get(handles.RemoveTPS_edit,'String'),' Images Finished.\n']);
    cd(WorkDir);
end

% Slice timing
if get(handles.SliceTiming_checkbox,'Value') == 1
    fprintf('Slice Timing Started.\n');
    ScanOrder = get(handles.SliceTiming_popupmenu,'Value');
    TR = str2double(get(handles.TR_edit,'String'));
    if get(handles.Parallel_checkbox,'Value') == 1
        parfor i=1:SubNum
            fprintf(['Slice Timing:',Sublist{i},' Started.\n']);
            spmjob = load([ProgramPath,filesep,'jobmats',filesep,'SliceTiming.mat']);
            cd([WorkDir,filesep,StartFolder,filesep,Sublist{i}]);
            Filename = dir('*.nii');
            Nii = nifti(Filename.name);
            FileList = cell(Nii.dat.dim(1,4),1);
            for j = 1:Nii.dat.dim(1,4)
                FileList{j} = [WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,Filename.name,',',num2str(j)];
            end
            spmjob.matlabbatch{1,1}.spm.temporal.st.scans{1} = FileList;
            nSlice = Nii.dat.dim(1,3);
            spmjob.matlabbatch{1,1}.spm.temporal.st.nslices = nSlice;
            switch ScanOrder
                case 1
                    spmjob.matlabbatch{1,1}.spm.temporal.st.so = [1:2:nSlice,2:2:nSlice];
                case 2
                    spmjob.matlabbatch{1,1}.spm.temporal.st.so = [2:2:nSlice,1:2:nSlice];
                case 3
                    if nSlice/2 == floor(nSlice/2)
                        spmjob.matlabbatch{1,1}.spm.temporal.st.so = [nSlice-1:-2:1,nSlice :-2:2];
                    else
                        spmjob.matlabbatch{1,1}.spm.temporal.st.so = [nSlice:-2:1,nSlice-1:-2:2];
                    end
                case 4
                    if nSlice/2 == floor(nSlice/2)
                        spmjob.matlabbatch{1,1}.spm.temporal.st.so = [nSlice:-2:2,nSlice-1:-2:1];
                    else
                        spmjob.matlabbatch{1,1}.spm.temporal.st.so = [nSlice-1:-2:2,nSlice:-2:1];
                    end
                case 5
                    spmjob.matlabbatch{1,1}.spm.temporal.st.so = 1:nSlice;
                case 6
                    spmjob.matlabbatch{1,1}.spm.temporal.st.so = nSlice:-1:1;
            end
            spmjob.matlabbatch{1,1}.spm.temporal.st.refslice = spmjob.matlabbatch{1,1}.spm.temporal.st.so(ceil(nSlice/2));
            spmjob.matlabbatch{1,1}.spm.temporal.st.tr = TR;
            spmjob.matlabbatch{1,1}.spm.temporal.st.ta = TR - (TR/nSlice);
            spm_jobman('initcfg');
            spm_jobman('run',spmjob.matlabbatch);
            mkdir([WorkDir,filesep,StartFolder,'A',filesep,Sublist{i}]);
            movefile([WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,'a*.nii'],[WorkDir,filesep,StartFolder,'A',filesep,Sublist{i}]);
            fprintf(['Slice Timing:',Sublist{i},' OK.\n']);
        end
    else
        for i=1:SubNum
            fprintf(['Slice Timing:',Sublist{i},' Started.\n']);
            spmjob = load([ProgramPath,filesep,'jobmats',filesep,'SliceTiming.mat']);
            cd([WorkDir,filesep,StartFolder,filesep,Sublist{i}]);
            Filename = dir('*.nii');
            Nii = nifti(Filename.name);
            FileList = cell(Nii.dat.dim(1,4),1);
            for j = 1:Nii.dat.dim(1,4)
                FileList{j} = [WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,Filename.name,',',num2str(j)];
            end
            spmjob.matlabbatch{1,1}.spm.temporal.st.scans{1} = FileList;
            nSlice = Nii.dat.dim(1,3);
            spmjob.matlabbatch{1,1}.spm.temporal.st.nslices = nSlice;
            switch ScanOrder
                case 1
                    spmjob.matlabbatch{1,1}.spm.temporal.st.so = [1:2:nSlice,2:2:nSlice];
                case 2
                    spmjob.matlabbatch{1,1}.spm.temporal.st.so = [2:2:nSlice,1:2:nSlice];
                case 3
                    if nSlice/2 == floor(nSlice/2)
                        spmjob.matlabbatch{1,1}.spm.temporal.st.so = [nSlice-1:-2:1,nSlice :-2:2];
                    else
                        spmjob.matlabbatch{1,1}.spm.temporal.st.so = [nSlice:-2:1,nSlice-1:-2:2];
                    end
                case 4
                    if nSlice/2 == floor(nSlice/2)
                        spmjob.matlabbatch{1,1}.spm.temporal.st.so = [nSlice:-2:2,nSlice-1:-2:1];
                    else
                        spmjob.matlabbatch{1,1}.spm.temporal.st.so = [nSlice-1:-2:2,nSlice:-2:1];
                    end
                case 5
                    spmjob.matlabbatch{1,1}.spm.temporal.st.so = 1:nSlice;
                case 6
                    spmjob.matlabbatch{1,1}.spm.temporal.st.so = nSlice:-1:1;
            end
            spmjob.matlabbatch{1,1}.spm.temporal.st.refslice = spmjob.matlabbatch{1,1}.spm.temporal.st.so(ceil(nSlice/2));
            spmjob.matlabbatch{1,1}.spm.temporal.st.tr = TR;
            spmjob.matlabbatch{1,1}.spm.temporal.st.ta = TR - (TR/nSlice);
            spm_jobman('initcfg');
            spm_jobman('run',spmjob.matlabbatch);
            mkdir([WorkDir,filesep,StartFolder,'A',filesep,Sublist{i}]);
            movefile([WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,'a*.nii'],[WorkDir,filesep,StartFolder,'A',filesep,Sublist{i}]);
            fprintf(['Slice Timing:',Sublist{i},' OK.\n']);
        end
    end
    fprintf('Slice Timing Finished.\n');
    StartFolder = [StartFolder,'A'];
    cd(WorkDir);
end

% Realign
if get(handles.Realign_checkbox,'Value') == 1
    fprintf('Realign Started.\n');
    mkdir([WorkDir,filesep,'HeadMotionParameter']);
    MaxHM = zeros(SubNum,6);
    if get(handles.Parallel_checkbox,'Value') == 1
        parfor i = 1:SubNum
            fprintf(['Realign:',Sublist{i},' Started.\n']);
            spmjob = load([ProgramPath,filesep,'jobmats',filesep,'Realign.mat']);
            cd([WorkDir,filesep,StartFolder,filesep,Sublist{i}]);
            Filename = dir('*.nii');
            Nii = nifti(Filename.name);
            FileList = cell(Nii.dat.dim(1,4),1);
            for j = 1:Nii.dat.dim(1,4)
                FileList{j} = [WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,Filename.name,',',num2str(j)];
            end
            spmjob.matlabbatch{1,1}.spm.spatial.realign.estwrite.data{1} = FileList;
            spm_jobman('initcfg');
            spm_jobman('run',spmjob.matlabbatch);
            mkdir([WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i}]);
            mkdir([WorkDir,filesep,StartFolder,'R',filesep,Sublist{i}])
            movefile([WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,'mean*'],[WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i}]);
            rpname = dir([WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,'rp*.txt']);
            headmotion = load(rpname.name);
            movefile(rpname.name,[WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i}]);
            movefile([WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,'r*.nii'],[WorkDir,filesep,StartFolder,'R',filesep,Sublist{i}]);
            headmotionMX = max(abs(headmotion));
            headmotionMX(4:6) = headmotionMX(4:6)*180/pi;
            MaxHM(i,:) = headmotionMX;
            fprintf(['Realign:',Sublist{i},' OK.\n']);
        end
    else
        for i = 1:SubNum
            fprintf(['Realign:',Sublist{i},' Started.\n']);
            spmjob = load([ProgramPath,filesep,'jobmats',filesep,'Realign.mat']);
            cd([WorkDir,filesep,StartFolder,filesep,Sublist{i}]);
            Filename = dir('*.nii');
            Nii = nifti(Filename.name);
            FileList = cell(Nii.dat.dim(1,4),1);
            for j = 1:Nii.dat.dim(1,4)
                FileList{j} = [WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,Filename.name,',',num2str(j)];
            end
            spmjob.matlabbatch{1,1}.spm.spatial.realign.estwrite.data{1} = FileList;
            spm_jobman('initcfg');
            spm_jobman('run',spmjob.matlabbatch);
            mkdir([WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i}]);
            mkdir([WorkDir,filesep,StartFolder,'R',filesep,Sublist{i}])
            movefile([WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,'mean*'],[WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i}]);
            rpname = dir([WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,'rp*.txt']);
            headmotion = load(rpname.name);
            movefile(rpname.name,[WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i}]);
            movefile([WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,'r*.nii'],[WorkDir,filesep,StartFolder,'R',filesep,Sublist{i}]);
            headmotionMX = max(abs(headmotion));
            headmotionMX(4:6) = headmotionMX(4:6)*180/pi;
            MaxHM(i,:) = headmotionMX;
            fprintf(['Realign:',Sublist{i},' OK.\n']);
        end
    end
    ExcludeSub_Text=[];
    for ExcludingCriteria = 3:-0.5:0.5
        BigestHM = max(MaxHM,[],2);
        ind = find(BigestHM>ExcludingCriteria);
        TempText = '';
        for ExcludingSub = 1:length(ind)
            TempText=sprintf('%s%s\n',TempText,Sublist{ind(ExcludingSub)});
        end
        ExcludeSub_Text=sprintf('%s\nExcluding Criteria: %2.1fmm and %2.1f degree in max head motion\n%s\n\n\n',ExcludeSub_Text,ExcludingCriteria,ExcludingCriteria,TempText);
    end
    fid = fopen([WorkDir,filesep,'HeadMotionParameter',filesep,'SubExcludingMaxHM.txt'],'at+');
    fprintf(fid,'%s',ExcludeSub_Text);
    fclose(fid);
    fprintf('Realign Started.\n');
    StartFolder = [StartFolder,'R'];
    cd(WorkDir);
end

% Normalization
if get(handles.Normalize_checkbox,'Value') == 1
    fprintf('Normalization Started.\n');
    SPMPath = spm('Dir');
    VoxSize = eval(['[',get(handles.Normalize_edit,'String'),']']);
    SPMversion = spm('Ver');
    switch get(handles.Normalize_popupmenu,'Value')
        case 1
            % Using EPI template
            if get(handles.Parallel_checkbox,'Value') == 1
                parfor i = 1:SubNum
                    fprintf(['Normalize:',Sublist{i},' Started.\n']);
                    cd([WorkDir,filesep,StartFolder,filesep,Sublist{i}]);
                    Filename = dir('*.nii');
                    Nii = nifti(Filename.name);
                    FileList = cell(Nii.dat.dim(1,4),1);
                    for j = 1:Nii.dat.dim(1,4)
                        FileList{j} = [WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,Filename.name,',',num2str(j)];
                    end
                    DirMean = dir([WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i},filesep,'mean*.nii']);
                    MeanFilename = [WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i},filesep,DirMean(1).name];
                    spmjob = load([ProgramPath,filesep,'jobmats',filesep,'Normalize.mat']);
                    spmjob.matlabbatch{1,1}.spm.spatial.normalise.estwrite.subj(1,1).source = {MeanFilename};
                    spmjob.matlabbatch{1,1}.spm.spatial.normalise.estwrite.subj(1,1).resample = [FileList;{MeanFilename}];
                    spmjob.matlabbatch{1,1}.spm.spatial.normalise.estwrite.eoptions.template = {[spm('Dir'),filesep,'templates',filesep,'EPI.nii,1']};
                    spmjob.matlabbatch{1,1}.spm.spatial.normalise.estwrite.roptions.bb = [-90,-126,-72;90,90,108];
                    spmjob.matlabbatch{1,1}.spm.spatial.normalise.estwrite.roptions.vox = VoxSize;
                    if str2double(SPMversion(4:end)) > 8
                        oldnorm = spmjob.matlabbatch{1,1}.spm.spatial.normalise;
                        oldnorm.estwrite.eoptions.template = {[SPMPath,filesep,'toolbox',filesep,'OldNorm',filesep,'EPI.nii,1']};
                        spmjob = [];
                        spmjob.matlabbatch{1,1}.spm.tools.oldnorm = oldnorm;
                    end
                    spm_jobman('initcfg');
                    spm_jobman('run',spmjob.matlabbatch);
                    mkdir([WorkDir,filesep,StartFolder,'W',filesep,Sublist{i}])
                    movefile([WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,'w*.nii'],[WorkDir,filesep,StartFolder,'W',filesep,Sublist{i}]);
                    fprintf(['Normalization:',Sublist{i},' OK.\n']);
                end
            else
                for i = 1:SubNum
                    fprintf(['Normalize:',Sublist{i},' Started.\n']);
                    cd([WorkDir,filesep,StartFolder,filesep,Sublist{i}]);
                    Filename = dir('*.nii');
                    Nii = nifti(Filename.name);
                    FileList = cell(Nii.dat.dim(1,4),1);
                    for j = 1:Nii.dat.dim(1,4)
                        FileList{j} = [WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,Filename.name,',',num2str(j)];
                    end
                    DirMean = dir([WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i},filesep,'mean*.nii']);
                    MeanFilename = [WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i},filesep,DirMean(1).name];
                    spmjob = load([ProgramPath,filesep,'jobmats',filesep,'Normalize.mat']);
                    spmjob.matlabbatch{1,1}.spm.spatial.normalise.estwrite.subj(1,1).source = {MeanFilename};
                    spmjob.matlabbatch{1,1}.spm.spatial.normalise.estwrite.subj(1,1).resample = [FileList;{MeanFilename}];
                    spmjob.matlabbatch{1,1}.spm.spatial.normalise.estwrite.eoptions.template = {[spm('Dir'),filesep,'templates',filesep,'EPI.nii,1']};
                    spmjob.matlabbatch{1,1}.spm.spatial.normalise.estwrite.roptions.bb = [-90,-126,-72;90,90,108];
                    spmjob.matlabbatch{1,1}.spm.spatial.normalise.estwrite.roptions.vox = VoxSize;
                    if str2double(SPMversion(4:end)) > 8
                        oldnorm = spmjob.matlabbatch{1,1}.spm.spatial.normalise;
                        oldnorm.estwrite.eoptions.template = {[SPMPath,filesep,'toolbox',filesep,'OldNorm',filesep,'EPI.nii,1']};
                        spmjob = [];
                        spmjob.matlabbatch{1,1}.spm.tools.oldnorm = oldnorm;
                    end
                    spm_jobman('initcfg');
                    spm_jobman('run',spmjob.matlabbatch);
                    mkdir([WorkDir,filesep,StartFolder,'W',filesep,Sublist{i}])
                    movefile([WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,'w*.nii'],[WorkDir,filesep,StartFolder,'W',filesep,Sublist{i}]);
                    fprintf(['Normalization:',Sublist{i},' OK.\n']);
                end
            end
        case 2
            % Using T1 segment
            if get(handles.Parallel_checkbox,'Value') == 1
                parfor i = 1:SubNum
                    % Coregister T1 to Fun
                    fprintf(['Coregester T1 to Fun:',Sublist{i},' Started.\n']);
                    mkdir([WorkDir,filesep,'T1ImgCoreg',filesep,Sublist{i}]);
                    copyfile([WorkDir,filesep,'T1Img',filesep,Sublist{i},filesep,'c*.nii'],[WorkDir,filesep,'T1ImgCoreg',filesep,Sublist{i}]);
                    DirMean = dir([WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i},filesep,'mean*.nii']);
                    spmjob = load([ProgramPath,filesep,'jobmats',filesep,'Coregister.mat']);
                    spmjob.matlabbatch{1,1}.spm.spatial.coreg.estimate.ref = {[WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i},filesep,DirMean.name,',1']};
                    SourceFile = dir([WorkDir,filesep,'T1ImgCoreg',filesep,Sublist{i},filesep,'*.nii']);
                    spmjob.matlabbatch{1,1}.spm.spatial.coreg.estimate.source={[WorkDir,filesep,'T1ImgCoreg',filesep,Sublist{i},filesep,SourceFile.name]};
                    fprintf(['Coregister:',Sublist{i},' OK.\n']);
                    spm_jobman('initcfg');
                    spm_jobman('run',spmjob.matlabbatch);
                    
                    % Segment
                    fprintf(['Segment:',Sublist{i},' Started.\n']);
                    mkdir([WorkDir,filesep,'T1ImgSegment',filesep,Sublist{i}]);
                    copyfile([WorkDir,filesep,'T1ImgCoreg',filesep,Sublist{i},filesep,'c*.nii'],[WorkDir,filesep,'T1ImgSegment',filesep,Sublist{i}]);
                    if str2double(SPMversion(4:end)) > 8
                        spmjob = load([ProgramPath,filesep,'jobmats',filesep,'Segment_12.mat']);
                        SourceFile = dir([WorkDir,filesep,'T1ImgSegment',filesep,Sublist{i},filesep,'*.nii']);
                        spmjob.matlabbatch{1,1}.spm.spatial.preproc.channel.vols = {[WorkDir,filesep,'T1ImgSegment',filesep,Sublist{i},filesep,SourceFile.name]};
                        for j = 1:6
                            spmjob.matlabbatch{1,1}.spm.spatial.preproc.tissue(j).tpm = {[SPMPath,filesep,'tpm',filesep,'TPM.nii,',num2str(j)]};
                        end
                    else
                        spmjob = load([ProgramPath,filesep,'jobmats',filesep,'Segment.mat']);
                        spmjob.matlabbatch{1,1}.spm.spatial.preproc.opts.tpm = {[SPMPath,filesep,'tpm',filesep,'grey.nii'];[SPMPath,filesep,'tpm',filesep,'white.nii'];[SPMPath,filesep,'tpm',filesep,'csf.nii']};
                        SourceFile = dir([WorkDir,filesep,'T1ImgSegment',filesep,Sublist{i},filesep,'*.nii']);
                        spmjob.matlabbatch{1,1}.spm.spatial.preproc.data={[WorkDir,filesep,'T1ImgSegment',filesep,Sublist{i},filesep,SourceFile.name]};
                        spmjob.matlabbatch{1,1}.spm.spatial.preproc.opts.regtype = 'mni';
                        oldseg = spmjob.matlabbatch{1,1}.spm.spatial.preproc;
                        oldseg.opts.tpm = {[SPMPath,filesep,'toolbox',filesep,'OldSeg',filesep,'grey.nii'];[SPMPath,filesep,'toolbox',filesep,'OldSeg',filesep,'white.nii'];[SPMPath,filesep,'toolbox',filesep,'OldSeg',filesep,'csf.nii']};
                        spmjob = [];
                        spmjob.matlabbatch{1,1}.spm.tools.oldseg = oldseg;
                    end
                    fprintf(['Segment:',Sublist{i},' OK.\n']);
                    spm_jobman('initcfg');
                    spm_jobman('run',spmjob.matlabbatch);
                    
                    % Normalize functional images
                    fprintf(['Normalize-Write:',Sublist{i},' Started.\n']);
                    Filename = dir([WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,'*.nii']);
                    Nii = nifti([WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,Filename.name]);
                    FileList = cell(Nii.dat.dim(1,4),1);
                    for j = 1:Nii.dat.dim(1,4)
                        FileList{j} = [WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,Filename.name,',',num2str(j)];
                    end
                    DirMean = dir([WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i},filesep,'mean*.nii']);
                    MeanFilename = [WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i},filesep,DirMean(1).name];
                    
                    if str2double(SPMversion(4:end)) > 8
                        spmjob = load([ProgramPath,filesep,'jobmats',filesep,'Normalize_Write_12.mat']);
                        MatFilename = dir([WorkDir,filesep,'T1ImgSegment',filesep,Sublist{i},filesep,'y_*.nii']);
                        spmjob.matlabbatch{1,1}.spm.spatial.normalise.write.subj.def = {[WorkDir,filesep,'T1ImgSegment',filesep,Sublist{i},filesep,MatFilename.name]};
                        spmjob.matlabbatch{1,1}.spm.spatial.normalise.write.subj.resample = [FileList;{MeanFilename}];
                    else
                        spmjob = load([ProgramPath,filesep,'jobmats',filesep,'Normalize_Write.mat']);
                        MatFilename = dir([WorkDir,filesep,'T1ImgSegment',filesep,Sublist{i},filesep,'*seg_sn.mat']);
                        spmjob.matlabbatch{1,1}.spm.spatial.normalise.write.subj(1,1).matname = {[WorkDir,filesep,'T1ImgSegment',filesep,Sublist{i},filesep,MatFilename.name]};
                        spmjob.matlabbatch{1,1}.spm.spatial.normalise.write.subj(1,1).resample = [FileList;{MeanFilename}];
                        spmjob.matlabbatch{1,1}.spm.spatial.normalise.write.roptions.bb = [-90,-126,-72;90,90,108];
                        spmjob.matlabbatch{1,1}.spm.spatial.normalise.write.roptions.vox = VoxSize;
                        oldnorm = spmjob.matlabbatch{1,1}.spm.spatial.normalise;
                        spmjob=[];
                        spmjob.matlabbatch{1,1}.spm.tools.oldnorm = oldnorm;
                    end
                    spm_jobman('initcfg');
                    spm_jobman('run',spmjob.matlabbatch);
                    mkdir([WorkDir,filesep,StartFolder,'W',filesep,Sublist{i}])
                    movefile([WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,'w*.nii'],[WorkDir,filesep,StartFolder,'W',filesep,Sublist{i}]);
                    fprintf(['Normalize-Write:',Sublist{i},' OK.\n']);
                end
            else
                for i = 1:SubNum
                    % Coregister T1 to Fun
                    fprintf(['Coregester T1 to Fun:',Sublist{i},' Started.\n']);
                    mkdir([WorkDir,filesep,'T1ImgCoreg',filesep,Sublist{i}]);
                    copyfile([WorkDir,filesep,'T1Img',filesep,Sublist{i},filesep,'c*.nii'],[WorkDir,filesep,'T1ImgCoreg',filesep,Sublist{i}]);
                    DirMean = dir([WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i},filesep,'mean*.nii']);
                    spmjob = load([ProgramPath,filesep,'jobmats',filesep,'Coregister.mat']);
                    spmjob.matlabbatch{1,1}.spm.spatial.coreg.estimate.ref = {[WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i},filesep,DirMean.name,',1']};
                    SourceFile = dir([WorkDir,filesep,'T1ImgCoreg',filesep,Sublist{i},filesep,'*.nii']);
                    spmjob.matlabbatch{1,1}.spm.spatial.coreg.estimate.source={[WorkDir,filesep,'T1ImgCoreg',filesep,Sublist{i},filesep,SourceFile.name]};
                    fprintf(['Coregister:',Sublist{i},' OK.\n']);
                    spm_jobman('initcfg');
                    spm_jobman('run',spmjob.matlabbatch);
                    
                    % Segment
                     fprintf(['Segment:',Sublist{i},' Started.\n']);
                    mkdir([WorkDir,filesep,'T1ImgSegment',filesep,Sublist{i}]);
                    copyfile([WorkDir,filesep,'T1ImgCoreg',filesep,Sublist{i},filesep,'c*.nii'],[WorkDir,filesep,'T1ImgSegment',filesep,Sublist{i}]);                    
                    if str2double(SPMversion(4:end)) > 8
                        spmjob = load([ProgramPath,filesep,'jobmats',filesep,'Segment_12.mat']);
                        SourceFile = dir([WorkDir,filesep,'T1ImgSegment',filesep,Sublist{i},filesep,'*.nii']);
                        spmjob.matlabbatch{1,1}.spm.spatial.preproc.channel.vols = {[WorkDir,filesep,'T1ImgSegment',filesep,Sublist{i},filesep,SourceFile.name]};
                        for j = 1:6
                            spmjob.matlabbatch{1,1}.spm.spatial.preproc.tissue(j).tpm = {[SPMPath,filesep,'tpm',filesep,'TPM.nii,',num2str(j)]};
                        end
                    else
                        spmjob = load([ProgramPath,filesep,'jobmats',filesep,'Segment.mat']);
                        spmjob.matlabbatch{1,1}.spm.spatial.preproc.opts.tpm = {[SPMPath,filesep,'tpm',filesep,'grey.nii'];[SPMPath,filesep,'tpm',filesep,'white.nii'];[SPMPath,filesep,'tpm',filesep,'csf.nii']};
                        SourceFile = dir([WorkDir,filesep,'T1ImgSegment',filesep,Sublist{i},filesep,'*.nii']);
                        spmjob.matlabbatch{1,1}.spm.spatial.preproc.data={[WorkDir,filesep,'T1ImgSegment',filesep,Sublist{i},filesep,SourceFile.name]};
                        spmjob.matlabbatch{1,1}.spm.spatial.preproc.opts.regtype = 'mni';
                        oldseg = spmjob.matlabbatch{1,1}.spm.spatial.preproc;
                        oldseg.opts.tpm = {[SPMPath,filesep,'toolbox',filesep,'OldSeg',filesep,'grey.nii'];[SPMPath,filesep,'toolbox',filesep,'OldSeg',filesep,'white.nii'];[SPMPath,filesep,'toolbox',filesep,'OldSeg',filesep,'csf.nii']};
                        spmjob = [];
                        spmjob.matlabbatch{1,1}.spm.tools.oldseg = oldseg;
                    end
                    fprintf(['Segment:',Sublist{i},' OK.\n']);
                    spm_jobman('initcfg');
                    spm_jobman('run',spmjob.matlabbatch);
                    
                    % Normalize functional images
                    fprintf(['Normalize-Write:',Sublist{i},' Started.\n']);
                    Filename = dir([WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,'*.nii']);
                    Nii = nifti([WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,Filename.name]);
                    FileList = cell(Nii.dat.dim(1,4),1);
                    for j = 1:Nii.dat.dim(1,4)
                        FileList{j} = [WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,Filename.name,',',num2str(j)];
                    end
                    DirMean = dir([WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i},filesep,'mean*.nii']);
                    MeanFilename = [WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i},filesep,DirMean(1).name];
                    
                    if str2double(SPMversion(4:end)) > 8
                        spmjob = load([ProgramPath,filesep,'jobmats',filesep,'Normalize_Write_12.mat']);
                        MatFilename = dir([WorkDir,filesep,'T1ImgSegment',filesep,Sublist{i},filesep,'y_*.nii']);
                        spmjob.matlabbatch{1,1}.spm.spatial.normalise.write.subj.def = {[WorkDir,filesep,'T1ImgSegment',filesep,Sublist{i},filesep,MatFilename.name]};
                        spmjob.matlabbatch{1,1}.spm.spatial.normalise.write.subj.resample = [FileList;{MeanFilename}];
                    else
                        spmjob = load([ProgramPath,filesep,'jobmats',filesep,'Normalize_Write.mat']);
                        MatFilename = dir([WorkDir,filesep,'T1ImgSegment',filesep,Sublist{i},filesep,'*seg_sn.mat']);
                        spmjob.matlabbatch{1,1}.spm.spatial.normalise.write.subj(1,1).matname = {[WorkDir,filesep,'T1ImgSegment',filesep,Sublist{i},filesep,MatFilename.name]};
                        spmjob.matlabbatch{1,1}.spm.spatial.normalise.write.subj(1,1).resample = [FileList;{MeanFilename}];
                        spmjob.matlabbatch{1,1}.spm.spatial.normalise.write.roptions.bb = [-90,-126,-72;90,90,108];
                        spmjob.matlabbatch{1,1}.spm.spatial.normalise.write.roptions.vox = VoxSize;
                        oldnorm = spmjob.matlabbatch{1,1}.spm.spatial.normalise;
                        spmjob=[];
                        spmjob.matlabbatch{1,1}.spm.tools.oldnorm = oldnorm;
                    end
                    spm_jobman('initcfg');
                    spm_jobman('run',spmjob.matlabbatch);
                    mkdir([WorkDir,filesep,StartFolder,'W',filesep,Sublist{i}])
                    movefile([WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,'w*.nii'],[WorkDir,filesep,StartFolder,'W',filesep,Sublist{i}]);
                    fprintf(['Normalize-Write:',Sublist{i},' OK.\n']);
                end
            end
    end
    
    
    % Generate Normalize Check Pic
    fprintf('Generating Picture for Normalization.\n');
    mkdir([WorkDir,filesep,'PicForChkNormalize']);
    cd([WorkDir,filesep,'PicForChkNormalize']);
    
    Ch2Filename = [ProgramPath,filesep,'..',filesep,'templates',filesep,'ch2.nii'];
    ch2_hdr = spm_vol(Ch2Filename);
    ch2_vol = spm_read_vols(ch2_hdr);
    
    c_view_u = ch2_vol(:,ceil(end/2),:);
    c_view_u = squeeze(c_view_u)';
    c_view_u = c_view_u(end:-1:1,:);
    c_view_u = c_view_u./max(c_view_u(:))*255;
    c_view_u = imresize(c_view_u,217/181);
    
    s_view_u = ch2_vol(ceil(end/2),:,:);
    s_view_u = squeeze(s_view_u)';
    s_view_u = s_view_u(end:-1:1,:);
    s_view_u = s_view_u./max(s_view_u(:))*255;
    s_view_u = imresize(s_view_u,217/181);
    
    a_view_u = ch2_vol(:,:,ceil(end/2));
    a_view_u = squeeze(a_view_u)';
    a_view_u = a_view_u(end:-1:1,:);
    a_view_u = a_view_u./max(a_view_u(:))*255;
    
    underlay = uint8([c_view_u,s_view_u,a_view_u]);
    underlay = repmat(underlay,[1,1,3]);
    
    if get(handles.Parallel_checkbox,'Value') == 1
        parfor i = 1:SubNum
            DirMean = dir([WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i},filesep,'wmean*.nii']);
            MeanFilename = [WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i},filesep,DirMean(1).name];
            
            mean_hdr = spm_vol(MeanFilename);
            mean_vol = spm_read_vols(mean_hdr);
            
            c_view_o = mean_vol(:,ceil(end/2),:);
            c_view_o = squeeze(c_view_o)';
            c_view_o = c_view_o(end:-1:1,:);
            c_view_o = c_view_o./max(c_view_o(:))*255;
            c_view_o = imresize(c_view_o,size(c_view_u));
            
            s_view_o = mean_vol(ceil(end/2),:,:);
            s_view_o = squeeze(s_view_o)';
            s_view_o = s_view_o(end:-1:1,:);
            s_view_o = s_view_o./max(s_view_o(:))*255;
            s_view_o = imresize(s_view_o,size(s_view_u));
            
            a_view_o = mean_vol(:,:,ceil(end/2));
            a_view_o = squeeze(a_view_o)';
            a_view_o = a_view_o(end:-1:1,:);
            a_view_o = a_view_o./max(a_view_o(:))*255;
            a_view_o = imresize(a_view_o,size(a_view_u));
            
            overlay = uint8([c_view_o,s_view_o,a_view_o]);
            overlay = repmat(overlay,[1,1,3]);
            overlay(:,:,2:3) = overlay(:,:,2:3)/2;
            outputimg = imresize(imadd(underlay./2,overlay./2),2);
            imwrite(outputimg,[Sublist{i},'.tif']);
        end
    else
        for i = 1:SubNum
            DirMean = dir([WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i},filesep,'wmean*.nii']);
            MeanFilename = [WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i},filesep,DirMean(1).name];
            
            mean_hdr = spm_vol(MeanFilename);
            mean_vol = spm_read_vols(mean_hdr);
            
            c_view_o = mean_vol(:,ceil(end/2),:);
            c_view_o = squeeze(c_view_o)';
            c_view_o = c_view_o(end:-1:1,:);
            c_view_o = c_view_o./max(c_view_o(:))*255;
            c_view_o = imresize(c_view_o,size(c_view_u));
            
            s_view_o = mean_vol(ceil(end/2),:,:);
            s_view_o = squeeze(s_view_o)';
            s_view_o = s_view_o(end:-1:1,:);
            s_view_o = s_view_o./max(s_view_o(:))*255;
            s_view_o = imresize(s_view_o,size(s_view_u));
            
            a_view_o = mean_vol(:,:,ceil(end/2));
            a_view_o = squeeze(a_view_o)';
            a_view_o = a_view_o(end:-1:1,:);
            a_view_o = a_view_o./max(a_view_o(:))*255;
            a_view_o = imresize(a_view_o,size(a_view_u));
            
            overlay = uint8([c_view_o,s_view_o,a_view_o]);
            overlay = repmat(overlay,[1,1,3]);
            overlay(:,:,2:3) = overlay(:,:,2:3)/2;
            outputimg = imresize(imadd(underlay./2,overlay./2),2);
            imwrite(outputimg,[Sublist{i},'.tif']);
        end
    end
    StartFolder = [StartFolder,'W'];
    fprintf('Normalization Finished.\n');
    cd(WorkDir);
end




% Smooth
if get(handles.Smooth_checkbox,'Value') == 1
    fprintf('Smoothing Started.\n');
    FWHM = eval(['[',get(handles.Smooth_edit,'String'),']']);
    if get(handles.Parallel_checkbox,'Value') == 1
        parfor i = 1:SubNum
            fprintf(['Smoothing:',Sublist{i},' Started.\n']);
            spmjob = load([ProgramPath,filesep,'jobmats',filesep,'Smooth.mat']);
            Filename = dir([WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,'*.nii']);
            Nii = nifti([WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,Filename.name]);
            FileList = cell(Nii.dat.dim(1,4),1);
            for j = 1:Nii.dat.dim(1,4)
                FileList{j} = [WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,Filename.name,',',num2str(j)];
            end
            spmjob.matlabbatch{1,1}.spm.spatial.smooth.data = FileList;
            spmjob.matlabbatch{1,1}.spm.spatial.smooth.fwhm = FWHM;
            spm_jobman('initcfg');
            spm_jobman('run',spmjob.matlabbatch);
            mkdir([WorkDir,filesep,StartFolder,'S',filesep,Sublist{i}])
            movefile([WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,'s*.nii'],[WorkDir,filesep,StartFolder,'S',filesep,Sublist{i}]);
            fprintf(['Smoothing:',Sublist{i},' OK.\n']);
        end
    else
        for i = 1:SubNum
            fprintf(['Smoothing:',Sublist{i},' Started.\n']);
            spmjob = load([ProgramPath,filesep,'jobmats',filesep,'Smooth.mat']);
            Filename = dir([WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,'*.nii']);
            Nii = nifti([WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,Filename.name]);
            FileList = cell(Nii.dat.dim(1,4),1);
            for j = 1:Nii.dat.dim(1,4)
                FileList{j} = [WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,Filename.name,',',num2str(j)];
            end
            spmjob.matlabbatch{1,1}.spm.spatial.smooth.data = FileList;
            spmjob.matlabbatch{1,1}.spm.spatial.smooth.fwhm = FWHM;
            spm_jobman('initcfg');
            spm_jobman('run',spmjob.matlabbatch);
            mkdir([WorkDir,filesep,StartFolder,'S',filesep,Sublist{i}])
            movefile([WorkDir,filesep,StartFolder,filesep,Sublist{i},filesep,'s*.nii'],[WorkDir,filesep,StartFolder,'S',filesep,Sublist{i}]);
            fprintf(['Smoothing:',Sublist{i},' OK.\n']);
        end
    end
    StartFolder = [StartFolder,'S'];
    cd(WorkDir);
    fprintf('Smoothing Finished.\n');
end

% Detrend
if get(handles.Detrend_checkbox,'Value') == 1
    fprintf('Detrend Started.\n');
    if get(handles.Parallel_checkbox,'Value') == 1
        parfor i = 1:SubNum
            fprintf(['Detrend:',Sublist{i},' Started.\n']);
            cd([WorkDir,filesep,StartFolder,filesep,Sublist{i}])
            Filename = dir('*.nii');
            Nii = nifti(Filename.name);
            tmpCourse = reshape(double(Nii.dat),[Nii.dat.dim(1,1) * Nii.dat.dim(1,2) * Nii.dat.dim(1,3),Nii.dat.dim(1,4)])';
            VoxelNum = size(tmpCourse,2);
            DetrendedCourse = zeros(size(tmpCourse));
            for j = 1:DataSegNum
                if j ~= DataSegNum
                    ind = ceil(VoxelNum/DataSegNum)*(j-1)+1:ceil(VoxelNum/DataSegNum)*j;
                else
                    ind = ceil(VoxelNum/DataSegNum)*(j-1)+1:VoxelNum;
                end
                DetrendedCourse(:,ind) = detrend(tmpCourse(:,ind)) + repmat(mean(tmpCourse(:,ind)),Nii.dat.dim(1,4),1);
            end
            DetrendedCourse = reshape(DetrendedCourse',[Nii.dat.dim(1,1),Nii.dat.dim(1,2),Nii.dat.dim(1,3),Nii.dat.dim(1,4)]);
            Nii.dat.fname = ['d',Filename.name];
            Nii.dat.dtype  = 16;
            create(Nii);
            Nii.dat(:,:,:,:) = DetrendedCourse;
            mkdir([WorkDir,filesep,StartFolder,'D',filesep,Sublist{i}])
            movefile('d*.nii',[WorkDir,filesep,StartFolder,'D',filesep,Sublist{i}]);
            fprintf(['Detrend:',Sublist{i},' OK.\n']);
        end
    else
        for i = 1:SubNum
            fprintf(['Detrend:',Sublist{i},' Started.\n']);
            cd([WorkDir,filesep,StartFolder,filesep,Sublist{i}])
            Filename = dir('*.nii');
            Nii = nifti(Filename.name);
            tmpCourse = reshape(double(Nii.dat),[Nii.dat.dim(1,1) * Nii.dat.dim(1,2) * Nii.dat.dim(1,3),Nii.dat.dim(1,4)])';
            VoxelNum = size(tmpCourse,2);
            DetrendedCourse = zeros(size(tmpCourse));
            for j = 1:DataSegNum
                if j ~= DataSegNum
                    ind = ceil(VoxelNum/DataSegNum)*(j-1)+1:ceil(VoxelNum/DataSegNum)*j;
                else
                    ind = ceil(VoxelNum/DataSegNum)*(j-1)+1:VoxelNum;
                end
                DetrendedCourse(:,ind) = detrend(tmpCourse(:,ind)) + repmat(mean(tmpCourse(:,ind)),Nii.dat.dim(1,4),1);
            end
            DetrendedCourse = reshape(DetrendedCourse',[Nii.dat.dim(1,1),Nii.dat.dim(1,2),Nii.dat.dim(1,3),Nii.dat.dim(1,4)]);
            Nii.dat.fname = ['d',Filename.name];
            Nii.dat.dtype  = 16;
            create(Nii);
            Nii.dat(:,:,:,:) = DetrendedCourse;
            mkdir([WorkDir,filesep,StartFolder,'D',filesep,Sublist{i}])
            movefile('d*.nii',[WorkDir,filesep,StartFolder,'D',filesep,Sublist{i}]);
            fprintf(['Detrend:',Sublist{i},' OK.\n']);
        end
    end
    StartFolder = [StartFolder,'D'];
    cd(WorkDir);
    fprintf('Detrend Finished.\n');
end

% Regress Covariates
if get(handles.Regress_checkbox,'Value') == 1
    fprintf('Regressing Out Covariates Started.\n');
    isHM = get(handles.RegHM_checkbox,'Value');
    HMtype = get(handles.RegHM_popupmenu,'Value');
    isWM = get(handles.RegWM_checkbox,'Value');
    isCSF = get(handles.RegCSF_checkbox,'Value');
    isGS = get(handles.RegGS_checkbox,'Value');
    if get(handles.Parallel_checkbox,'Value') == 1
        parfor i = 1:SubNum
            fprintf(['Regressing Out Covariates:',Sublist{i},' Started.\n']);
            HMCov = [];
            WMCov = [];
            CSFCov = [];
            GSCov = [];
            if isHM == 1
                rpname = dir([WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i},filesep,'rp*.txt']);
                headmotion  = load([WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i},filesep,rpname.name]);
                switch HMtype
                    case 1
                        HMCov = headmotion;
                    case 2
                        HMCov = [headmotion, [zeros(1,size(headmotion,2));headmotion(1:end-1,:)], headmotion.^2, [zeros(1,size(headmotion,2));headmotion(1:end-1,:)].^2];
                end
            end
            cd([WorkDir,filesep,StartFolder,filesep,Sublist{i}])
            Filename = dir('*.nii');
            Nii = nifti(Filename.name);
            tmpCourse = reshape(double(Nii.dat),[Nii.dat.dim(1,1) * Nii.dat.dim(1,2) * Nii.dat.dim(1,3),Nii.dat.dim(1,4)])';
            if isWM == 1
                WM_hdr = spm_vol([ProgramPath,filesep,'..',filesep,'templates',filesep,'WhiteMask_09_61x73x61.nii']);
                WM_vol = spm_read_vols(WM_hdr);
                WMCov = mean(tmpCourse(:,reshape(WM_vol>0,1,[])),2);
            end
            if isCSF == 1
                CSF_hdr = spm_vol([ProgramPath,filesep,'..',filesep,'templates',filesep,'CsfMask_07_61x73x61.nii']);
                CSF_vol = spm_read_vols(CSF_hdr);
                CSFCov = mean(tmpCourse(:,reshape(CSF_vol>0,1,[])),2);
            end
            if isGS == 1
                GS_hdr = spm_vol([ProgramPath,filesep,'..',filesep,'templates',filesep,'BrainMask_05_61x73x61.nii']);
                GS_vol = spm_read_vols(GS_hdr);
                GSCov = mean(tmpCourse(:,reshape(GS_vol>0,1,[])),2);
            end
            AllCov = [ones(Nii.dat.dim(1,4),1),HMCov,WMCov,CSFCov,GSCov];
            
            VoxelNum = size(tmpCourse,2);
            RegressedCourse = zeros(size(tmpCourse));
            for j = 1:DataSegNum
                if j ~= DataSegNum
                    ind = ceil(VoxelNum/DataSegNum)*(j-1)+1:ceil(VoxelNum/DataSegNum)*j;
                else
                    ind = ceil(VoxelNum/DataSegNum)*(j-1)+1:VoxelNum;
                end
                RegressedCourse(:,ind) = RegressCov(tmpCourse(:,ind),AllCov);
            end
            RegressedCourse(isnan(RegressedCourse)) = 0;
            RegressedCourse = reshape(RegressedCourse',[Nii.dat.dim(1,1),Nii.dat.dim(1,2),Nii.dat.dim(1,3),Nii.dat.dim(1,4)]);
            Nii.dat.fname = ['c',Filename.name];
            Nii.dat.dtype  = 16;
            create(Nii);
            Nii.dat(:,:,:,:) = RegressedCourse;
            mkdir([WorkDir,filesep,StartFolder,'C',filesep,Sublist{i}])
            movefile('c*.nii',[WorkDir,filesep,StartFolder,'C',filesep,Sublist{i}]);
            fprintf(['Regressing Out Covariates:',Sublist{i},' OK.\n']);
        end
    else
        for i = 1:SubNum
            fprintf(['Regressing Out Covariates:',Sublist{i},' Started.\n']);
            HMCov = [];
            WMCov = [];
            CSFCov = [];
            GSCov = [];
            if isHM == 1
                rpname = dir([WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i},filesep,'rp*.txt']);
                headmotion  = load([WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i},filesep,rpname.name]);
                switch HMtype
                    case 1
                        HMCov = headmotion;
                    case 2
                        HMCov = [headmotion, [zeros(1,size(headmotion,2));headmotion(1:end-1,:)], headmotion.^2, [zeros(1,size(headmotion,2));headmotion(1:end-1,:)].^2];
                end
            end
            cd([WorkDir,filesep,StartFolder,filesep,Sublist{i}])
            Filename = dir('*.nii');
            Nii = nifti(Filename.name);
            tmpCourse = reshape(double(Nii.dat),[Nii.dat.dim(1,1) * Nii.dat.dim(1,2) * Nii.dat.dim(1,3),Nii.dat.dim(1,4)])';
            if isWM == 1
                WM_hdr = spm_vol([ProgramPath,filesep,'..',filesep,'templates',filesep,'WhiteMask_09_61x73x61.nii']);
                WM_vol = spm_read_vols(WM_hdr);
                WMCov = mean(tmpCourse(:,reshape(WM_vol>0,1,[])),2);
            end
            if isCSF == 1
                CSF_hdr = spm_vol([ProgramPath,filesep,'..',filesep,'templates',filesep,'CsfMask_07_61x73x61.nii']);
                CSF_vol = spm_read_vols(CSF_hdr);
                CSFCov = mean(tmpCourse(:,reshape(CSF_vol>0,1,[])),2);
            end
            if isGS == 1
                GS_hdr = spm_vol([ProgramPath,filesep,'..',filesep,'templates',filesep,'BrainMask_05_61x73x61.nii']);
                GS_vol = spm_read_vols(GS_hdr);
                GSCov = mean(tmpCourse(:,reshape(GS_vol>0,1,[])),2);
            end
            AllCov = [ones(Nii.dat.dim(1,4),1),HMCov,WMCov,CSFCov,GSCov];
            
            VoxelNum = size(tmpCourse,2);
            RegressedCourse = zeros(size(tmpCourse));
            for j = 1:DataSegNum
                if j ~= DataSegNum
                    ind = ceil(VoxelNum/DataSegNum)*(j-1)+1:ceil(VoxelNum/DataSegNum)*j;
                else
                    ind = ceil(VoxelNum/DataSegNum)*(j-1)+1:VoxelNum;
                end
                RegressedCourse(:,ind) = RegressCov(tmpCourse(:,ind),AllCov);
            end
            RegressedCourse(isnan(RegressedCourse)) = 0;
            RegressedCourse = reshape(RegressedCourse',[Nii.dat.dim(1,1),Nii.dat.dim(1,2),Nii.dat.dim(1,3),Nii.dat.dim(1,4)]);
            Nii.dat.fname = ['c',Filename.name];
            Nii.dat.dtype  = 16;
            create(Nii);
            Nii.dat(:,:,:,:) = RegressedCourse;
            mkdir([WorkDir,filesep,StartFolder,'C',filesep,Sublist{i}])
            movefile('c*.nii',[WorkDir,filesep,StartFolder,'C',filesep,Sublist{i}]);
            fprintf(['Regressing Out Covariates:',Sublist{i},' OK.\n']);
        end
    end
    StartFolder = [StartFolder,'C'];
    cd(WorkDir);
    fprintf('Regressing Out Covariates Finished.\n');
end

%Filter
if get(handles.Filter_checkbox,'Value') == 1
    fprintf('Filtering Started.\n');
    FilterLow = str2double(get(handles.FilterLow_edit,'String'));
    FilterHigh = str2double(get(handles.FilterHigh_edit,'String'));
    TR = str2double(get(handles.TR_edit,'String'));
    if get(handles.Parallel_checkbox,'Value') == 1
        parfor i = 1:SubNum
            fprintf(['Filtering:',Sublist{i},' Started.\n']);
            cd([WorkDir,filesep,StartFolder,filesep,Sublist{i}])
            Filename = dir('*.nii');
            Nii = nifti(Filename.name);
            tmpCourse = reshape(double(Nii.dat),[Nii.dat.dim(1,1) * Nii.dat.dim(1,2) * Nii.dat.dim(1,3),Nii.dat.dim(1,4)])';
            VoxelNum = size(tmpCourse,2);
            FilteredCourse = zeros(size(tmpCourse));
            for j = 1:DataSegNum
                if j ~= DataSegNum
                    ind = ceil(VoxelNum/DataSegNum)*(j-1)+1:ceil(VoxelNum/DataSegNum)*j;
                else
                    ind = ceil(VoxelNum/DataSegNum)*(j-1)+1:VoxelNum;
                end
                FilteredCourse(:,ind) = BandFilter(tmpCourse(:,ind),TR,[FilterLow,FilterHigh]);
            end
            FilteredCourse = reshape(FilteredCourse',[Nii.dat.dim(1,1),Nii.dat.dim(1,2),Nii.dat.dim(1,3),Nii.dat.dim(1,4)]);
            Nii.dat.fname = ['f',Filename.name];
            Nii.dat.dtype  = 16;
            create(Nii);
            Nii.dat(:,:,:,:) = FilteredCourse;
            mkdir([WorkDir,filesep,StartFolder,'F',filesep,Sublist{i}])
            movefile('f*.nii',[WorkDir,filesep,StartFolder,'F',filesep,Sublist{i}]);
            fprintf(['Filtering:',Sublist{i},' OK.\n']);
        end
    else
        for i = 1:SubNum
            fprintf(['Filtering:',Sublist{i},' Started.\n']);
            cd([WorkDir,filesep,StartFolder,filesep,Sublist{i}])
            Filename = dir('*.nii');
            Nii = nifti(Filename.name);
            tmpCourse = reshape(double(Nii.dat),[Nii.dat.dim(1,1) * Nii.dat.dim(1,2) * Nii.dat.dim(1,3),Nii.dat.dim(1,4)])';
            VoxelNum = size(tmpCourse,2);
            FilteredCourse = zeros(size(tmpCourse));
            for j = 1:DataSegNum
                if j ~= DataSegNum
                    ind = ceil(VoxelNum/DataSegNum)*(j-1)+1:ceil(VoxelNum/DataSegNum)*j;
                else
                    ind = ceil(VoxelNum/DataSegNum)*(j-1)+1:VoxelNum;
                end
                FilteredCourse(:,ind) = BandFilter(tmpCourse(:,ind),TR,[FilterLow,FilterHigh]);
            end
            FilteredCourse = reshape(FilteredCourse',[Nii.dat.dim(1,1),Nii.dat.dim(1,2),Nii.dat.dim(1,3),Nii.dat.dim(1,4)]);
            Nii.dat.fname = ['f',Filename.name];
            Nii.dat.dtype  = 16;
            create(Nii);
            Nii.dat(:,:,:,:) = FilteredCourse;
            mkdir([WorkDir,filesep,StartFolder,'F',filesep,Sublist{i}])
            movefile('f*.nii',[WorkDir,filesep,StartFolder,'F',filesep,Sublist{i}]);
            fprintf(['Filtering:',Sublist{i},' OK.\n']);
        end
    end
    StartFolder = [StartFolder,'F'];
    cd(WorkDir);
    fprintf('Filtering Finished.\n');
end

%Scrubbing
if get(handles.Scrubbing_checkbox,'Value') == 1
    fprintf('Scrubbing Started.\n');
    FDThr = str2double(get(handles.ScrubbingThr_edit,'String'));
    TpsBefore = str2double(get(handles.ScrubbingBefore_edit,'String'));
    TpsAfter = str2double(get(handles.ScrubbingAfter_edit,'String'));
    Scrub_meth = get(handles.Scrubbing_popupmenu,'Value');
    Scrub_percent = zeros(SubNum,1);
    Scrub_meanFD = zeros(SubNum,1);
    if get(handles.Parallel_checkbox,'Value') == 1
        parfor i = 1:SubNum
            fprintf(['Scrubbing:',Sublist{i},' Started.\n']);
            cd([WorkDir,filesep,StartFolder,filesep,Sublist{i}])
            Filename = dir('*.nii');
            Nii = nifti(Filename.name);
            tmpCourse = reshape(double(Nii.dat),[Nii.dat.dim(1,1) * Nii.dat.dim(1,2) * Nii.dat.dim(1,3),Nii.dat.dim(1,4)])';
            rpname = dir([WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i},filesep,'rp*.txt']);
            headmotion  = load([WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i},filesep,rpname.name]);
            FDPower = [zeros(1,6);diff(headmotion)];
            FDPower(:,4:6) = FDPower(:,4:6)*50;
            FDPower = sum(abs(FDPower),2);
            ind = find(FDPower>FDThr);
            Scrub_percent(i) = length(ind)/Nii.dat.dim(1,4);
            Scrub_meanFD(i) = mean(FDPower);
            FDIndex = ones(Nii.dat.dim(1,4),1);
            FDIndex(ind) = 0;
            for j = 1:TpsBefore
                indtmp = ind - j;
                indtmp(indtmp<1) = 1;
                FDIndex(indtmp) = 0;
            end
            for j = 1:TpsAfter
                indtmp = ind + j;
                indtmp(indtmp>Nii.dat.dim(1,4)) = Nii.dat.dim(1,4);
                FDIndex(indtmp) = 0;
            end
            ScrubbingCourse = tmpCourse(FDIndex>0,:);
            xi = 1:Nii.dat.dim(1,4);
            x = xi(FDIndex>0);
            VoxelNum = size(tmpCourse,2);
            if Scrub_meth ~= 1
                for j = 1:DataSegNum
                    if j ~= DataSegNum
                        ind = ceil(VoxelNum/DataSegNum)*(j-1)+1:ceil(VoxelNum/DataSegNum)*j;
                    else
                        ind = ceil(VoxelNum/DataSegNum)*(j-1)+1:VoxelNum;
                    end
                    tmpCourse(:,ind) = Scrubbing(ScrubbingCourse(:,ind),x,xi,Scrub_meth);
                end
                ScrubbingCourse = tmpCourse;
            end
            ntime = size(ScrubbingCourse,1);
            ScrubbingCourse = reshape(ScrubbingCourse',[Nii.dat.dim(1,1),Nii.dat.dim(1,2),Nii.dat.dim(1,3),ntime]);
            Nii.dat.fname = ['b',Filename.name];
            Nii.dat.dtype  = 16;
            Nii.dat.dim(1,4) = ntime;
            create(Nii);
            Nii.dat(:,:,:,:) = ScrubbingCourse;
            mkdir([WorkDir,filesep,StartFolder,'B',filesep,Sublist{i}])
            movefile('b*.nii',[WorkDir,filesep,StartFolder,'B',filesep,Sublist{i}]);
            fprintf(['Scrubbing:',Sublist{i},' OK.\n']);
        end
    else
        for i = 1:SubNum
            fprintf(['Scrubbing:',Sublist{i},' Started.\n']);
            cd([WorkDir,filesep,StartFolder,filesep,Sublist{i}])
            Filename = dir('*.nii');
            Nii = nifti(Filename.name);
            tmpCourse = reshape(double(Nii.dat),[Nii.dat.dim(1,1) * Nii.dat.dim(1,2) * Nii.dat.dim(1,3),Nii.dat.dim(1,4)])';
            rpname = dir([WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i},filesep,'rp*.txt']);
            headmotion  = load([WorkDir,filesep,'HeadMotionParameter',filesep,Sublist{i},filesep,rpname.name]);
            FDPower = [zeros(1,6);diff(headmotion)];
            FDPower(:,4:6) = FDPower(:,4:6)*50;
            FDPower = sum(abs(FDPower),2);
            ind = find(FDPower>FDThr);
            Scrub_percent(i) = length(ind)/Nii.dat.dim(1,4);
            Scrub_meanFD(i) = mean(FDPower);
            FDIndex = ones(Nii.dat.dim(1,4),1);
            FDIndex(ind) = 0;
            for j = 1:TpsBefore
                indtmp = ind - j;
                indtmp(indtmp<1) = 1;
                FDIndex(indtmp) = 0;
            end
            for j = 1:TpsAfter
                indtmp = ind + j;
                indtmp(indtmp>Nii.dat.dim(1,4)) = Nii.dat.dim(1,4);
                FDIndex(indtmp) = 0;
            end
            ScrubbingCourse = tmpCourse(FDIndex>0,:);
            xi = 1:Nii.dat.dim(1,4);
            x = xi(FDIndex>0);
            VoxelNum = size(tmpCourse,2);
            if Scrub_meth ~= 1
                for j = 1:DataSegNum
                    if j ~= DataSegNum
                        ind = ceil(VoxelNum/DataSegNum)*(j-1)+1:ceil(VoxelNum/DataSegNum)*j;
                    else
                        ind = ceil(VoxelNum/DataSegNum)*(j-1)+1:VoxelNum;
                    end
                    tmpCourse(:,ind) = Scrubbing(ScrubbingCourse(:,ind),x,xi,Scrub_meth);
                end
                ScrubbingCourse = tmpCourse;
            end
            ntime = size(ScrubbingCourse,1);
            ScrubbingCourse = reshape(ScrubbingCourse',[Nii.dat.dim(1,1),Nii.dat.dim(1,2),Nii.dat.dim(1,3),ntime]);
            Nii.dat.fname = ['b',Filename.name];
            Nii.dat.dtype  = 16;
            Nii.dat.dim(1,4) = ntime;
            create(Nii);
            Nii.dat(:,:,:,:) = ScrubbingCourse;
            mkdir([WorkDir,filesep,StartFolder,'B',filesep,Sublist{i}])
            movefile('b*.nii',[WorkDir,filesep,StartFolder,'B',filesep,Sublist{i}]);
            fprintf(['Scrubbing:',Sublist{i},' OK.\n']);
        end
    end
    ScrubPer_Text='ID   Scrub_Percentage   Mean_FD';
    for i = 1:SubNum
        ScrubPer_Text = sprintf('%s\n%s  %1.4f  %1.4f',ScrubPer_Text,Sublist{i},Scrub_percent(i),Scrub_meanFD(i));
    end
    fid = fopen([WorkDir,filesep,'HeadMotionParameter',filesep,'ScrubbingPercentage.txt'],'at+');
    fprintf(fid,'%s',ScrubPer_Text);
    fclose(fid);
    cd(WorkDir);
    fprintf('Scrubbing Finished.\n');
end


toc
set(hObject,'Enable','on');
set(hObject,'String','RUN');

function Data_Scrubbed = Scrubbing(Data,x,xi,meth)
switch meth
    case 2
        Data_Scrubbed = interp1(x,Data,xi,'nearest','extrap');
    case 3
        Data_Scrubbed = interp1(x,Data,xi,'linear','extrap');
    case 4
        Data_Scrubbed = interp1(x,Data,xi,'spline','extrap');
end


function Data_Filtered = BandFilter(Data, SamplePeriod, Band)
sampleFreq 	 = 1/SamplePeriod;
sampleLength = size(Data,1);
paddedLength = 2^nextpow2(sampleLength);
LowCutoff_HighPass = Band(1);
HighCutoff_LowPass = Band(2);
if (LowCutoff_HighPass >= sampleFreq/2)
    idxLowCutoff_HighPass = paddedLength/2 + 1;
else
    idxLowCutoff_HighPass = ceil(LowCutoff_HighPass * paddedLength * SamplePeriod + 1);
end

if (HighCutoff_LowPass>=sampleFreq/2)||(HighCutoff_LowPass==0)s
    idxHighCutoff_LowPass = paddedLength/2 + 1;
else
    idxHighCutoff_LowPass = fix(HighCutoff_LowPass * paddedLength * SamplePeriod + 1);
end
FrequencyMask = zeros(paddedLength,1);
FrequencyMask(idxLowCutoff_HighPass:idxHighCutoff_LowPass,1) = 1;
FrequencyMask(paddedLength-idxLowCutoff_HighPass+2:-1:paddedLength-idxHighCutoff_LowPass+2,1) = 1;
FrequencySetZero_Index = FrequencyMask==0;
Data = Data - repmat(mean(Data),size(Data,1),1);
Data = [Data;zeros(paddedLength -sampleLength,size(Data,2))];
Data = fft(Data);
Data(FrequencySetZero_Index,:) = 0;
Data = ifft(Data);
Data_Filtered = Data(1:sampleLength,:);

function r = RegressCov(y,X)
[n,ncolX] = size(X);
[Q,R,perm] = qr(X,0);
p = sum(abs(diag(R)) > max(n,ncolX)*eps(R(1)));
if p < ncolX,
    R = R(1:p,1:p);
    Q = Q(:,1:p);
    perm = perm(1:p);
end
b = zeros(ncolX,size(y,2));
b(perm,:) = R \ (Q'*y);
yhat = X * b;
r = y - yhat;

function dcm2nii(Filename, OutputDir)
if ispc
    filesep = '\';
else
    filesep = '/';
end

[~,~,ext] = fileparts(Filename);
if strcmpi(ext,'.nii')
    copyfile(Filename,OutputDir);
else
    ProgramPath = fileparts(which('SeeCAT_PrepfMRI.m'));
    OldDirTemp=pwd;
    cd([ProgramPath,filesep,'dcm2nii']);
    if ispc
        eval(['!dcm2nii.exe ','-b dcm2nii.ini',' -o ',OutputDir,' ',Filename]);
    elseif ismac
        eval(['!chmod +x dcm2nii_mac']);
        eval(['!./dcm2nii_mac ','-b dcm2nii.ini',' -o ',OutputDir,' ',Filename]);
    else
        eval('!chmod +x dcm2nii_linux');
        eval(['!./dcm2nii_linux ','-b dcm2nii.ini',' -o ',OutputDir,' ',Filename]);
    end
    cd(OldDirTemp);
end


% --- Executes on button press in Quit_pushbutton.
function Quit_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Quit_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.SeeCAT_PrefMRI_figure);


% --- Executes on button press in Detrend_checkbox.
function Detrend_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Detrend_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Detrend_checkbox


% --- Executes on button press in Regress_checkbox.
function Regress_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Regress_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Regress_checkbox
if get(hObject,'Value') == 1
    set(handles.RegHM_checkbox,'Enable','on');
    set(handles.RegWM_checkbox,'Enable','on');
    set(handles.RegCSF_checkbox,'Enable','on');
    set(handles.RegGS_checkbox,'Enable','on');
    if get(handles.RegHM_checkbox,'Value') == 1
        set(handles.RegHM_popupmenu,'Enable','on');
    else
        set(handles.RegHM_popupmenu,'Enable','off');
    end
else
    set(handles.RegHM_checkbox,'Enable','off');
    set(handles.RegHM_popupmenu,'Enable','off');
    set(handles.RegWM_checkbox,'Enable','off');
    set(handles.RegCSF_checkbox,'Enable','off');
    set(handles.RegGS_checkbox,'Enable','off');
end



% --- Executes on button press in RegHM_checkbox.
function RegHM_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to RegHM_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RegHM_checkbox
if get(hObject,'Value') == 1
    set(handles.RegHM_popupmenu,'Enable','on');
else
    set(handles.RegHM_popupmenu,'Enable','off');
end


% --- Executes on selection change in RegHM_popupmenu.
function RegHM_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to RegHM_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns RegHM_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from RegHM_popupmenu


% --- Executes during object creation, after setting all properties.
function RegHM_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RegHM_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RegWM_checkbox.
function RegWM_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to RegWM_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RegWM_checkbox


% --- Executes on button press in RegGS_checkbox.
function RegGS_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to RegGS_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RegGS_checkbox


% --- Executes on button press in RegCSF_checkbox.
function RegCSF_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to RegCSF_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RegCSF_checkbox


% --- Executes on button press in Filter_checkbox.
function Filter_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Filter_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Filter_checkbox
if get(hObject,'Value') == 1
    set(handles.FilterLow_edit,'Enable','on');
    set(handles.FilterHigh_edit,'Enable','on');
else
    set(handles.FilterLow_edit,'Enable','off');
    set(handles.FilterHigh_edit,'Enable','off');
end



function FilterLow_edit_Callback(hObject, eventdata, handles)
% hObject    handle to FilterLow_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FilterLow_edit as text
%        str2double(get(hObject,'String')) returns contents of FilterLow_edit as a double


% --- Executes during object creation, after setting all properties.
function FilterLow_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FilterLow_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FilterHigh_edit_Callback(hObject, eventdata, handles)
% hObject    handle to FilterHigh_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FilterHigh_edit as text
%        str2double(get(hObject,'String')) returns contents of FilterHigh_edit as a double


% --- Executes during object creation, after setting all properties.
function FilterHigh_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FilterHigh_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Realign_checkbox.
function Realign_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Realign_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Realign_checkbox


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in Normalize_checkbox.
function Normalize_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Normalize_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Normalize_checkbox
if get(hObject,'Value') == 1
    set(handles.Normalize_popupmenu,'Enable','on');
    set(handles.Normalize_edit,'Enable','on');
else
    set(handles.Normalize_popupmenu,'Enable','off');
    set(handles.Normalize_edit,'Enable','off');
end


% --- Executes on selection change in Normalize_popupmenu.
function Normalize_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to Normalize_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Normalize_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Normalize_popupmenu


% --- Executes during object creation, after setting all properties.
function Normalize_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Normalize_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Smooth_checkbox.
function Smooth_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Smooth_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Smooth_checkbox
if get(hObject,'Value') == 1
    set(handles.Smooth_edit,'Enable','on');
else
    set(handles.Smooth_edit,'Enable','off');
end


% --- Executes on button press in DICOMtoNIFTI_checkbox.
function DICOMtoNIFTI_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to DICOMtoNIFTI_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DICOMtoNIFTI_checkbox


% --- Executes on button press in RemoveTPS_checkbox.
function RemoveTPS_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveTPS_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RemoveTPS_checkbox
if get(hObject,'Value') == 1
    set(handles.RemoveTPS_edit,'Enable','on');
else
    set(handles.RemoveTPS_edit,'Enable','off');
end


% --- Executes on button press in SliceTiming_checkbox.
function SliceTiming_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to SliceTiming_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SliceTiming_checkbox
if get(hObject,'Value') == 1
    set(handles.SliceTiming_popupmenu,'Enable','on');
else
    set(handles.SliceTiming_popupmenu,'Enable','off');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Smooth_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Smooth_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Smooth_edit as text
%        str2double(get(hObject,'String')) returns contents of Smooth_edit as a double


% --- Executes during object creation, after setting all properties.
function Smooth_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Smooth_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Normalize_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Normalize_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Normalize_edit as text
%        str2double(get(hObject,'String')) returns contents of Normalize_edit as a double


% --- Executes during object creation, after setting all properties.
function Normalize_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Normalize_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ScrubbingThr_edit_Callback(hObject, eventdata, handles)
% hObject    handle to ScrubbingThr_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ScrubbingThr_edit as text
%        str2double(get(hObject,'String')) returns contents of ScrubbingThr_edit as a double


% --- Executes during object creation, after setting all properties.
function ScrubbingThr_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ScrubbingThr_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Scrubbing_checkbox.
function Scrubbing_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Scrubbing_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Scrubbing_checkbox
if get(hObject,'Value') == 1
    set(handles.ScrubbingThr_edit,'Enable','on');
    set(handles.Scrubbing_popupmenu,'Enable','on');
    set(handles.ScrubbingBefore_edit,'Enable','on');
    set(handles.ScrubbingAfter_edit,'Enable','on');
else
    set(handles.ScrubbingThr_edit,'Enable','off');
    set(handles.Scrubbing_popupmenu,'Enable','off');
    set(handles.ScrubbingBefore_edit,'Enable','off');
    set(handles.ScrubbingAfter_edit,'Enable','off');
end



function RemoveTPS_edit_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveTPS_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RemoveTPS_edit as text
%        str2double(get(hObject,'String')) returns contents of RemoveTPS_edit as a double


% --- Executes during object creation, after setting all properties.
function RemoveTPS_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RemoveTPS_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in SliceTiming_popupmenu.
function SliceTiming_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to SliceTiming_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SliceTiming_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SliceTiming_popupmenu


% --- Executes during object creation, after setting all properties.
function SliceTiming_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SliceTiming_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Parallel_checkbox.
function Parallel_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Parallel_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Parallel_checkbox



function ScrubbingBefore_edit_Callback(hObject, eventdata, handles)
% hObject    handle to ScrubbingBefore_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ScrubbingBefore_edit as text
%        str2double(get(hObject,'String')) returns contents of ScrubbingBefore_edit as a double


% --- Executes during object creation, after setting all properties.
function ScrubbingBefore_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ScrubbingBefore_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ScrubbingAfter_edit_Callback(hObject, eventdata, handles)
% hObject    handle to ScrubbingAfter_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ScrubbingAfter_edit as text
%        str2double(get(hObject,'String')) returns contents of ScrubbingAfter_edit as a double


% --- Executes during object creation, after setting all properties.
function ScrubbingAfter_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ScrubbingAfter_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Scrubbing_popupmenu.
function Scrubbing_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to Scrubbing_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Scrubbing_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Scrubbing_popupmenu


% --- Executes during object creation, after setting all properties.
function Scrubbing_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Scrubbing_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
