function varargout = SeeCAT_VoxDeg(varargin)
% SEECAT_VOXDEG MATLAB code for SeeCAT_VoxDeg.fig
%      SEECAT_VOXDEG, by itself, creates a new SEECAT_VOXDEG or raises the existing
%      singleton*.
%
%      H = SEECAT_VOXDEG returns the handle to a new SEECAT_VOXDEG or the handle to
%      the existing singleton*.
%
%      SEECAT_VOXDEG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEECAT_VOXDEG.M with the given input arguments.
%
%      SEECAT_VOXDEG('Property','Value',...) creates a new SEECAT_VOXDEG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SeeCAT_VoxDeg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SeeCAT_VoxDeg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SeeCAT_VoxDeg

% Last Modified by GUIDE v2.5 05-Jul-2017 12:25:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SeeCAT_VoxDeg_OpeningFcn, ...
    'gui_OutputFcn',  @SeeCAT_VoxDeg_OutputFcn, ...
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


% --- Executes just before SeeCAT_VoxDeg is made visible.
function SeeCAT_VoxDeg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SeeCAT_VoxDeg (see VARARGIN)

% Choose default command line output for SeeCAT_VoxDeg
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
movegui(hObject,'center');
% UIWAIT makes SeeCAT_VoxDeg wait for user response (see UIRESUME)
% uiwait(handles.figure1);
set(handles.WorkDir_edit,'String',pwd);
ProgramPath = fileparts(which('SeeCAT.m'));
if ispc
    filesep = '\';
else
    filesep = '/';
end
set(handles.Mask_edit,'String',[ProgramPath,filesep,'templates',filesep,'GreyMask_02_61x73x61.nii']);


% --- Outputs from this function are returned to the command line.
function varargout = SeeCAT_VoxDeg_OutputFcn(hObject, eventdata, handles)
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
MaskName = get(handles.Mask_edit,'String');


Log_Text = datestr(now);
Log_Text = sprintf('%s\n%s\n',Log_Text,['Working directory: ',WorkDir]);
Log_Text = sprintf('%s%s\n',Log_Text,['Start folder: ',StartFolder]);
Log_Text = sprintf('%s%s\n',Log_Text,['Calculation mask: ',MaskName]);
Log_Text = sprintf('%s%s\n',Log_Text,['Threshold: ',get(handles.Threshold_edit,'String')]);
if get(handles.Local_checkbox,'Value') == 1;
    Log_Text = sprintf('%s%s\n',Log_Text,['Exclude local connections: <',get(handles.Local_edit,'String'),' mm']);
end
Log_Text = sprintf('%s%s',Log_Text,'Output type: ');
if get(handles.Positive_checkbox,'Value') == 1;
    Log_Text = sprintf('%s%s',Log_Text,'Positive ');
end
if get(handles.Absolute_checkbox,'Value') == 1;
    Log_Text = sprintf('%s%s',Log_Text,'Absolute ');
end
if get(handles.Negative_checkbox,'Value') == 1;
    Log_Text = sprintf('%s%s',Log_Text,'Negative ');
end
Log_Text = sprintf('%s%s\n',Log_Text,'');
Log_Text = sprintf('%s%s',Log_Text,'Output metric: ');
if get(handles.Degree_checkbox,'Value') == 1;
    Log_Text = sprintf('%s%s',Log_Text,'Binarized ');
end
if get(handles.Strength_checkbox,'Value') == 1;
    Log_Text = sprintf('%s%s',Log_Text,'Weighted ');
end
Log_Text = sprintf('%s%s\n',Log_Text,'');
if get(handles.Distance_checkbox,'Value') == 1;
    Log_Text = sprintf('%s%s\n',Log_Text,['Distance: ',get(handles.Distance_edit,'String'),' mm']);
end
if get(handles.Normalize_checkbox,'Value') == 1;
    Log_Text = sprintf('%s%s\n',Log_Text,'Normalized to Z-score ');
end
if get(handles.Smooth_checkbox,'Value') == 1
    Log_Text = sprintf('%s%s\n',Log_Text,['Smooth: FWHM = [',get(handles.Smooth_edit,'String'),']']);
end

fid = fopen([WorkDir,filesep,'Log_VoxDeg_',datestr(now,'yyyymmddHHMMSS'),'.txt'],'at+');
fprintf(fid,'%s',Log_Text);
fclose(fid);

fprintf('Calculating Voxel-based Degree Started.\n');

hdr_mask = spm_vol(MaskName);
[vol_mask,XYZ] = spm_read_vols(hdr_mask);


Sublist = dir([WorkDir,filesep,StartFolder]);
if strcmpi(Sublist(3).name,'.DS_Store')
    Sublist(1:3) = [];
else
    Sublist(1:2) = [];
end
mask_ind = reshape(vol_mask>0,1,[]);
mask_ind2 = find(mask_ind);
mkdir([WorkDir,filesep,'Results',filesep,'VoxelDeg_',StartFolder]);
XYZ = XYZ(:,mask_ind);
CalType = [0,0,0];
if get(handles.Positive_checkbox,'Value') == 1
    CalType(1) = 1;
end
if get(handles.Absolute_checkbox,'Value') == 1
    CalType(2) = 1;
end
if get(handles.Negative_checkbox,'Value') == 1
    CalType(3) = 1;
end
nCalType = sum(CalType);

CalMatri = [0,0];
if get(handles.Degree_checkbox,'Value') == 1
    CalMatri(1) = 1;
end
if get(handles.Strength_checkbox,'Value') == 1
    CalMatri(2) = 1;
end
nCalMatri = sum(CalMatri);
if get(handles.Distance_checkbox,'Value') == 1
    CalDistanceSeg = [0,eval(['[',get(handles.Distance_edit,'String'),']']),inf];
else
    CalDistanceSeg = [0,inf];
end
nCalDistance = length(CalDistanceSeg) - 1;
if get(handles.Normalize_checkbox,'Value') == 1
    CalNormal = 1;
else
    CalNormal = 0;
end
if get(handles.Smooth_checkbox,'Value') == 1
    CalSmooth = 1;
    CalSmoothFWHM = eval(['[',get(handles.Smooth_edit,'String'),']']);
else
    CalSmooth = 0;
end

ii = num2cell(repmat((1:nCalType)',[1,nCalMatri,nCalDistance]));
jj = num2cell(repmat((1:nCalMatri),[nCalType,1,nCalDistance]));
k(1,1,:) = 1:nCalDistance;
kk = num2cell(repmat(k,[nCalType,nCalMatri,1]));

CalTypeCell{1} = CalType;
CalTypeCell = repmat(CalTypeCell,[nCalType,nCalMatri,nCalDistance]);
CalMatriCell{1} = CalMatri;
CalMatriCell = repmat(CalMatriCell,[nCalType,nCalMatri,nCalDistance]);
CalDistanceSegCell{1} = CalDistanceSeg;
CalDistanceSegCell = repmat(CalDistanceSegCell,[nCalType,nCalMatri,nCalDistance]);
name_VoxDeg = cellfun(@genname,ii,jj,kk,CalTypeCell,CalMatriCell,CalDistanceSegCell,'UniformOutput',false);
if get(handles.Local_checkbox,'Value') == 1
    exLocal = 1;
    exLocaldis = str2double(get(handles.Local_edit,'String'));
else
    exLocal = 0;
    exLocaldis = 0;
end
degthresh = str2double(get(handles.Threshold_edit,'String'));
hdr{1} = hdr_mask;
hdr = repmat(hdr,[size(name_VoxDeg,1),size(name_VoxDeg,2),size(name_VoxDeg,3)]);
if get(handles.Parallel_checkbox,'Value') == 1
    parfor i = 1:length(Sublist)
        if ~isdir([WorkDir,filesep,StartFolder,filesep,Sublist(i).name])
            continue;
        end
        fprintf(['Calculating Voxel-based Degree:',Sublist(i).name,' Started.\n']);
        vol_VoxDeg = cell(1);
        vol_VoxDeg{1} = zeros(hdr_mask.dim);
        vol_VoxDeg = repmat(vol_VoxDeg,[nCalType,nCalMatri,nCalDistance]);
        cd([WorkDir,filesep,StartFolder,filesep,Sublist(i).name]);
        Filename = dir('*.nii');
        Nii = nifti(Filename.name);
        volCourse = reshape(double(Nii.dat),[Nii.dat.dim(1,1) * Nii.dat.dim(1,2) * Nii.dat.dim(1,3),Nii.dat.dim(1,4)])';
        maskCourse = volCourse(:,mask_ind);
        maskCourse = maskCourse - repmat(mean(maskCourse),[Nii.dat.dim(1,4),1]);
        maskCourse = maskCourse./repmat(std(maskCourse,0,1),[Nii.dat.dim(1,4),1]);
        for j = 1:size(maskCourse,2)
            D = pdist2(XYZ(:,j)',XYZ');
            r = maskCourse(:,j)' * maskCourse / (Nii.dat.dim(1,4) - 1);
            for k = 1:nCalType
                for l = 1:nCalMatri
                    for m = 1:nCalDistance
                        r_tmp = r;
                        ind = find(CalType,k);
                        ind = ind(end);
                        switch ind
                            case 1
                                r_tmp(r_tmp<0) = 0;
                            case 2
                                r_tmp = abs(r_tmp);
                            case 3
                                r_tmp = -r_tmp;
                                r_tmp(r_tmp<0) = 0;
                        end
                        if exLocal == 1
                            r_tmp(D<exLocaldis) = 0;
                        end
                        r_tmp(r_tmp<degthresh) = 0;
                        
                        ind = find(CalMatri,l);
                        ind = ind(end);
                        if ind == 1
                            r_tmp(r_tmp>0) = 1;
                        end
                        r_tmp(D<=CalDistanceSeg(m)|D>CalDistanceSeg(m+1)) = 0;
                        r_tmp(isnan(r_tmp)) = 0;
                        vol_VoxDeg{k,l,m}(mask_ind2(j)) = sum(r_tmp);
                    end
                end
            end
        end
        cd([WorkDir,filesep,'Results',filesep,'VoxelDeg_',StartFolder]);
        Subname = cell(1);
        Subname{1} = Sublist(i).name;
        Subname = repmat(Subname,[size(name_VoxDeg,1),size(name_VoxDeg,2),size(name_VoxDeg,3)]);
        cellfun(@writedegvol,hdr,name_VoxDeg,Subname,vol_VoxDeg);
        
        if CalNormal == 1
            zvol_VoxDeg = vol_VoxDeg;
            for j = 1:numel(vol_VoxDeg)
                zvol_VoxDeg{j}(mask_ind) = (vol_VoxDeg{j}(mask_ind) - mean(vol_VoxDeg{j}(mask_ind)))./std(vol_VoxDeg{j}(mask_ind));
            end
            cellfun(@zwritedegvol,hdr,name_VoxDeg,Subname,zvol_VoxDeg);
        end
        fprintf(['Calculating Voxel-based Degree:',Sublist(i).name,' OK.\n']);
    end
else
    for i = 1:length(Sublist)
        if ~isdir([WorkDir,filesep,StartFolder,filesep,Sublist(i).name])
            continue;
        end
        fprintf(['Calculating Voxel-based Degree:',Sublist(i).name,' Started.\n']);
        vol_VoxDeg = cell(1);
        vol_VoxDeg{1} = zeros(hdr_mask.dim);
        vol_VoxDeg = repmat(vol_VoxDeg,[nCalType,nCalMatri,nCalDistance]);
        cd([WorkDir,filesep,StartFolder,filesep,Sublist(i).name]);
        Filename = dir('*.nii');
        Nii = nifti(Filename.name);
        volCourse = reshape(double(Nii.dat),[Nii.dat.dim(1,1) * Nii.dat.dim(1,2) * Nii.dat.dim(1,3),Nii.dat.dim(1,4)])';
        maskCourse = volCourse(:,mask_ind);
        maskCourse = maskCourse - repmat(mean(maskCourse),[Nii.dat.dim(1,4),1]);
        maskCourse = maskCourse./repmat(std(maskCourse,0,1),[Nii.dat.dim(1,4),1]);
        for j = 1:size(maskCourse,2)
            D = pdist2(XYZ(:,j)',XYZ');
            r = maskCourse(:,j)' * maskCourse / (Nii.dat.dim(1,4) - 1);
            for k = 1:nCalType
                for l = 1:nCalMatri
                    for m = 1:nCalDistance
                        r_tmp = r;
                        ind = find(CalType,k);
                        ind = ind(end);
                        switch ind
                            case 1
                                r_tmp(r_tmp<0) = 0;
                            case 2
                                r_tmp = abs(r_tmp);
                            case 3
                                r_tmp = -r_tmp;
                                r_tmp(r_tmp<0) = 0;
                        end
                        if exLocal == 1
                            r_tmp(D<exLocaldis) = 0;
                        end
                        r_tmp(r_tmp<degthresh) = 0;
                        
                        ind = find(CalMatri,l);
                        ind = ind(end);
                        if ind == 1
                            r_tmp(r_tmp>0) = 1;
                        end
                        r_tmp(D<=CalDistanceSeg(m)|D>CalDistanceSeg(m+1)) = 0;
                        r_tmp(isnan(r_tmp)) = 0;
                        vol_VoxDeg{k,l,m}(mask_ind2(j)) = sum(r_tmp);
                    end
                end
            end
        end
        cd([WorkDir,filesep,'Results',filesep,'VoxelDeg_',StartFolder]);
        Subname = cell(1);
        Subname{1} = Sublist(i).name;
        Subname = repmat(Subname,[size(name_VoxDeg,1),size(name_VoxDeg,2),size(name_VoxDeg,3)]);
        cellfun(@writedegvol,hdr,name_VoxDeg,Subname,vol_VoxDeg);
        
        if CalNormal == 1
            zvol_VoxDeg = vol_VoxDeg;
            for j = 1:numel(vol_VoxDeg)
                zvol_VoxDeg{j}(mask_ind) = (vol_VoxDeg{j}(mask_ind) - mean(vol_VoxDeg{j}(mask_ind)))./std(vol_VoxDeg{j}(mask_ind));
            end
            cellfun(@zwritedegvol,hdr,name_VoxDeg,Subname,zvol_VoxDeg);
        end
        fprintf(['Calculating Voxel-based Degree:',Sublist(i).name,' OK.\n']);
    end
end
cd([WorkDir,filesep,'Results',filesep,'VoxelDeg_',StartFolder]);
if CalSmooth == 1
    spmjob = load([fileparts(which('SeeCAT_PrepfMRI.m')),filesep,'jobmats',filesep,'Smooth.mat']);
    list = dir('*.nii');
    FileList = cell(length(list),1);
    for i = 1:length(list)
        FileList{i} = [WorkDir,filesep,'Results',filesep,'VoxelDeg_',StartFolder,filesep,list(i).name,',1'];
    end
    spmjob.matlabbatch{1,1}.spm.spatial.smooth.data = FileList;
    spmjob.matlabbatch{1,1}.spm.spatial.smooth.fwhm = CalSmoothFWHM;
    spm_jobman('initcfg');
    spm_jobman('run',spmjob.matlabbatch);
end

if CalSmooth == 1 && CalNormal == 1
    for i = 1:numel(name_VoxDeg)
        movefile(['sz',name_VoxDeg{i},'*.nii'],[WorkDir,filesep,'Results',filesep,'VoxelDeg_',StartFolder,filesep,'sz',name_VoxDeg{i}(1:end-1)]);
    end
end
if CalSmooth == 1
    for i = 1:numel(name_VoxDeg)
        movefile(['s',name_VoxDeg{i},'*.nii'],[WorkDir,filesep,'Results',filesep,'VoxelDeg_',StartFolder,filesep,'s',name_VoxDeg{i}(1:end-1)]);
    end
end
if CalNormal == 1
    for i = 1:numel(name_VoxDeg)
        movefile(['z',name_VoxDeg{i},'*.nii'],[WorkDir,filesep,'Results',filesep,'VoxelDeg_',StartFolder,filesep,'z',name_VoxDeg{i}(1:end-1)]);
    end
end
for i = 1:numel(name_VoxDeg)
    movefile([name_VoxDeg{i},'*.nii'],[WorkDir,filesep,'Results',filesep,'VoxelDeg_',StartFolder,filesep,name_VoxDeg{i}(1:end-1)]);
end
cd(WorkDir);
fprintf('Calculating Voxel-based Degree Finished.\n');
toc
set(hObject,'Enable','on');
set(hObject,'String','RUN');

function writedegvol(hdr,name,subname,vol)
hdr.fname = [name,subname,'.nii'];
hdr.dt(1) = 16;
spm_write_vol(hdr,vol);

function zwritedegvol(hdr,name,subname,vol)
hdr.fname = ['z',name,subname,'.nii'];
hdr.dt(1) = 16;
spm_write_vol(hdr,vol);

function name = genname(i,j,k,CalType,CalMatri,CalDistanceSeg)
ind = find(CalType,i);
ind = ind(end);
switch ind
    case 1
        name_type = 'Pos_';
    case 2
        name_type = 'Abs_';
    case 3
        name_type = 'Neg_';
end
ind = find(CalMatri,j);
ind = ind(end);
switch ind
    case 1
        name_metri = 'Bin_';
    case 2
        name_metri = 'Wei_';
end
startDis = ['00',num2str(CalDistanceSeg(k))];
startDis = startDis(end-2:end);
endDis = ['00',num2str(CalDistanceSeg(k+1))];
endDis = endDis(end-2:end);
name = [name_type,name_metri,startDis,'-',endDis,'_'];





% --- Executes on button press in Quit_pushbutton.
function Quit_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Quit_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);


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



function Threshold_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Threshold_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Threshold_edit as text
%        str2double(get(hObject,'String')) returns contents of Threshold_edit as a double


% --- Executes during object creation, after setting all properties.
function Threshold_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Threshold_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Local_checkbox.
function Local_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Local_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Local_checkbox
if get(hObject,'Value') == 1
    set(handles.Local_edit,'Enable','on');
else
    set(handles.Local_edit,'Enable','off');
end



function Local_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Local_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Local_edit as text
%        str2double(get(hObject,'String')) returns contents of Local_edit as a double


% --- Executes during object creation, after setting all properties.
function Local_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Local_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Degree_checkbox.
function Degree_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Degree_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Degree_checkbox



% --- Executes on button press in Strength_checkbox.
function Strength_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Strength_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Strength_checkbox


% --- Executes on button press in Positive_checkbox.
function Positive_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Positive_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Positive_checkbox


% --- Executes on button press in Absolute_checkbox.
function Absolute_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Absolute_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Absolute_checkbox


% --- Executes on button press in Negative_checkbox.
function Negative_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Negative_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Negative_checkbox


% --- Executes on button press in Distance_checkbox.
function Distance_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Distance_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Distance_checkbox
if get(hObject,'Value') == 1
    set(handles.Distance_edit,'Enable','on');
else
    set(handles.Distance_edit,'Enable','off');
end



function Distance_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Distance_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Distance_edit as text
%        str2double(get(hObject,'String')) returns contents of Distance_edit as a double


% --- Executes during object creation, after setting all properties.
function Distance_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Distance_edit (see GCBO)
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


% --- Executes on button press in Normalize_checkbox.
function Normalize_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Normalize_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Normalize_checkbox


% --- Executes on button press in Parallel_checkbox.
function Parallel_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Parallel_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Parallel_checkbox
