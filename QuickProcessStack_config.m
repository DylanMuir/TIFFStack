function varargout = QuickProcessStack_config(varargin)
% QUICKPROCESSSTACK_CONFIG MATLAB code for QuickProcessStack_config.fig
%      QUICKPROCESSSTACK_CONFIG, by itself, creates a new QUICKPROCESSSTACK_CONFIG or raises the existing
%      singleton*.
%
%      H = QUICKPROCESSSTACK_CONFIG returns the handle to a new QUICKPROCESSSTACK_CONFIG or the handle to
%      the existing singleton*.
%
%      QUICKPROCESSSTACK_CONFIG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in QUICKPROCESSSTACK_CONFIG.M with the given input arguments.
%
%      QUICKPROCESSSTACK_CONFIG('Property','Value',...) creates a new QUICKPROCESSSTACK_CONFIG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before QuickProcessStack_config_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to QuickProcessStack_config_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help QuickProcessStack_config

% Last Modified by GUIDE v2.5 28-Oct-2011 14:08:40


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @QuickProcessStack_config_OpeningFcn, ...
                   'gui_OutputFcn',  @QuickProcessStack_config_OutputFcn, ...
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


% --- Executes just before QuickProcessStack_config is made visible.
function QuickProcessStack_config_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*DEFNU,*INUSL,*INUSD>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to QuickProcessStack_config (see VARARGIN)

% Choose default command line output for QuickProcessStack_config
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(handles.strExtractionFun, 'ColumnFormat', {'char'}, 'ColumnEditable', true, ...
   'ColumnWidth', {'auto'}, 'ColumnName', {'Extraction function'}, 'ColumnWidth', {225}, ...
   'Data', {''});
set(handles.strUseStimIDs, 'ColumnFormat', {'char'}, 'ColumnEditable', true, ...
   'ColumnWidth', {'auto'}, 'ColumnName', {'StimIDs to use'}, 'ColumnWidth', {115}, ...
   'Data', {'2:9'});


% UIWAIT makes QuickProcessStack_config wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = QuickProcessStack_config_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pbSaveConfig.
function pbSaveConfig_Callback(hObject, eventdata, handles)
% hObject    handle to pbSaveConfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% - Extract configuration data
bAlign = get(handles.bAlign, 'Value'); %#ok<NASGU>
bAssignBlack = get(handles.bAssignBlack, 'Value'); %#ok<NASGU>
bAssignBlank = get(handles.bAssignBlank, 'Value'); %#ok<NASGU>
tFrameDuration = get(handles.tFrameDuration, 'Data'); %#ok<NASGU>
nBlankStimID = get(handles.nBlankStimID, 'Data'); %#ok<NASGU>
strUseStimIDs = get(handles.strUseStimIDs, 'Data'); %#ok<NASGU>
strExtractionFun = get(handles.strExtractionFun, 'Data'); %#ok<NASGU>
mfStimTimeData = get(handles.StimConfigTable, 'Data'); %#ok<NASGU>

% - Prompt for a file to save data, and save configuration
uisave({'bAlign', 'bAssignBlack', 'bAssignBlank', 'tFrameDuration', 'nBlankStimID', 'strUseStimIDs', 'strExtractionFun', 'mfStimTimeData'}, ...
       'QPS_config.mat');


% --- Executes on button press in pbLoadConfig.
function pbLoadConfig_Callback(hObject, eventdata, handles)
% hObject    handle to pbLoadConfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% - Prompt for a file to load configuration from
strConfigFilename = uigetfile('QPS_config.mat', 'Load a configuration file', 'QPS_config.mat');

% - Load data
sConfigData = load(strConfigFilename);

% - Check for the correct data
if (~all(isfield(sConfigData, {'bAlign', 'bAssignBlack', 'bAssignBlank', 'tFrameDuration', 'nBlankStimID', 'strUseStimIDs', 'strExtractionFun', 'mfStimTimeData'})))
   errordlg('Invalid saved configuration file', 'Invalid saved config', 'modal');
   
else
   % - Assign the data to the dialog
   set(handles.bAlign, 'Value', sConfigData.bAlign);
   set(handles.bAssignBlack, 'Value', sConfigData.bAssignBlack);
   set(handles.bAssignBlank, 'Value', sConfigData.bAssignBlank);
   set(handles.tFrameDuration, 'Data', sConfigData.tFrameDuration);
   set(handles.nBlankStimID, 'Data', sConfigData.nBlankStimID);
   set(handles.strUseStimIDs, 'Data', sConfigData.strUseStimIDs);
   set(handles.strExtractionFun, 'Data', sConfigData.strExtractionFun);
   set(handles.StimConfigTable, 'Data', sConfigData.mfStimTimeData);
end

% --- Executes on button press in pbOK.
function pbOK_Callback(hObject, eventdata, handles)
% hObject    handle to pbOK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% -- Global configuration variable
global QPS_sConfig;

% -- Assign configuration data to QPS_sConfig
QPS_sConfig.bAlign = get(handles.bAlign, 'Value');
QPS_sConfig.bAssignBlack = get(handles.bAssignBlack, 'Value');
QPS_sConfig.bAssignBlank = get(handles.bAssignBlank, 'Value');
QPS_sConfig.tFrameDuration = get(handles.tFrameDuration, 'Data');
QPS_sConfig.nBlankStimID = get(handles.nBlankStimID, 'Data');

% - Try to get stim IDs to use
try
   strUseStimIDs = get(handles.strUseStimIDs, 'Data');
   if (isempty(strUseStimIDs{1}))
      QPS_sConfig.vnUseStimIDs = [];
   else
      QPS_sConfig.vnUseStimIDs = eval(strUseStimIDs{1});
   end
catch %#ok<CTCH>
   errordlg('Invalid entry for stim IDs to use', 'Invalid configuration', 'modal');
   return;
end

% - Try to get an extraction function
try
   strExtractionFun = get(handles.strExtractionFun, 'Data');
   if (isempty(strExtractionFun{1}))
      QPS_sConfig.fhExtract = [];
   else
      QPS_sConfig.fhExtract = eval(strExtractionFun{1});
   end
catch %#ok<CTCH>
   errordlg('Invalid entry for extraction function', 'Invalid configuration', 'modal');
   return;
end

% - Extract info from stimulus time configuration table
cmfStimTimeData = get(handles.StimConfigTable, 'Data');

% - Identify filled cells, numeric cells
mbFilledCell = ~cellfun(@isempty, cmfStimTimeData);
mbNumericCell = cellfun(@isnumeric, cmfStimTimeData);
mbIsNan = false(size(cmfStimTimeData));
mbIsNan(mbFilledCell) = cellfun(@isnan, cmfStimTimeData(mbFilledCell));
vbIgnoreRow = all(mbIsNan, 2);

% - Check for filled rows
vbPartialRow = (any(mbFilledCell, 2) & ~all(mbFilledCell, 2));

if (any(vbPartialRow))
   errordlg('Rows defining stimuli in the stimulus time configuration table must have all entries filled.', ...
      'Invalid configuration', 'modal');
   return;
end

% - Check for numeric data
if (any(mbFilledCell(:) & ~mbNumericCell(:)))
   errordlg('The stimulus time configuration dialog can only contain numbers and NaN.', ...
      'Invalid configuration', 'modal');
   return;
end

% - Extract filled rows in the table
mbTakeCells = logical(mbFilledCell .* repmat(~vbIgnoreRow, 1, size(mbFilledCell, 2)));
mfStimTimeData = reshape(cell2mat(cmfStimTimeData(mbTakeCells)), [], 7);

QPS_sConfig.nNumStimuli = size(mfStimTimeData, 1);
QPS_sConfig.vtStimulusDurations = mfStimTimeData(:, 1);
QPS_sConfig.mtBlankTimes = mfStimTimeData(:, 2:3);
QPS_sConfig.mtStimulusUseTimes = mfStimTimeData(:, 4:5);
QPS_sConfig.mtStimulusLabelTimes = mfStimTimeData(:, 6:7);

% - Close dialog
delete(QuickProcessStack_config);
drawnow;



% --- Executes on button press in pbCancel.
function pbCancel_Callback(hObject, eventdata, handles)
% hObject    handle to pbCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(QuickProcessStack_config);

% --- Executes on button press in pbLoadRegions.
function pbLoadRegions_Callback(hObject, eventdata, handles)
% hObject    handle to pbLoadRegions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.strROIFilename, 'String', uigetfile('*.zip'));


% --- Executes during object creation, after setting all properties.
function strROIFilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to strROIFilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in bAlign.
function bAlign_Callback(hObject, eventdata, handles)
% hObject    handle to bAlign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bAlign


% --- Executes on button press in bAssignBlack.
function bAssignBlack_Callback(hObject, eventdata, handles)
% hObject    handle to bAssignBlack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bAssignBlack


% --- Executes on button press in bAssignBlank.
function bAssignBlank_Callback(hObject, eventdata, handles)
% hObject    handle to bAssignBlank (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bAssignBlank

% --- END of QuickProcessStack_config.m ---

