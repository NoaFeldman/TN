function varargout = PasteDataDialog(varargin)
% PASTEDATADIALOG MATLAB code for PasteDataDialog.fig
%      PASTEDATADIALOG, by itself, creates a new PASTEDATADIALOG or raises the existing
%      singleton*.
%
%      H = PASTEDATADIALOG returns the handle to a new PASTEDATADIALOG or the handle to
%      the existing singleton*.
%
%      PASTEDATADIALOG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PASTEDATADIALOG.M with the given input arguments.
%
%      PASTEDATADIALOG('Property','Value',...) creates a new PASTEDATADIALOG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PasteDataDialog_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PasteDataDialog_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PasteDataDialog

% Last Modified by GUIDE v2.5 28-Oct-2012 23:20:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @PasteDataDialog_OpeningFcn, ...
    'gui_OutputFcn',  @PasteDataDialog_OutputFcn, ...
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


% --- Executes just before PasteDataDialog is made visible.
function PasteDataDialog_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PasteDataDialog (see VARARGIN)

% Choose default command line output for PasteDataDialog
handles.output = hObject;
handles.Data = PasteDataDialogDataClass();
   
% Update handles structure
guidata(hObject, handles);

dataText = clipboard('paste');

try
    [data, status] = str2num(dataText);
    
    columnNames = {};
    
    if (~status)
        % Maybe the first line is headers?
        dataLines = textscan(dataText, '%s', 'Delimiter', '\n');
        dataLines = dataLines{1};
        
        if (numel(dataLines) > 1)
            headerLine = dataLines{1};
            
            if (~isempty(headerLine))
                columnNames = textscan(headerLine, '%s', 'Delimiter', ' \t');
                columnNames = strrep(columnNames{1}, char(26), 'd'); % Replace the delta symbol
            end
            
            for i = 2:numel(dataLines)
                dataLineNumbers = str2num(dataLines{i});
                data(end+1, :) = dataLineNumbers;
            end
            
            status = 1;
        end
    end
    
    if (status)
        set(handles.dataTable, 'Data', data);
        
        columnNumbers = [1:size(data, 2)]';
        columnNumberStrings = cellstr(num2str(columnNumbers));
        
        if (numel(columnNames) < numel(columnNumberStrings))
            columnNames(numel(columnNames)+1:numel(columnNumberStrings)) = columnNumberStrings(numel(columnNames)+1:numel(columnNumberStrings));
        else
            set(handles.dataTable, 'ColumnName', columnNames);
        end
        
        set(handles.xColPopup, 'String', columnNames);
        set(handles.yColPopup, 'String', columnNames);
        set(handles.dxColPopup, 'String', columnNames);
        set(handles.dyColPopup, 'String', columnNames);
        
        set(handles.xColPopup, 'Value', min(1, numel(columnNames)));
        set(handles.dxColPopup, 'Value', min(2, numel(columnNames)));
        set(handles.yColPopup, 'Value', min(3, numel(columnNames)));
        set(handles.dyColPopup, 'Value', min(4, numel(columnNames)));
    end
    
    if (isempty(data))
        dlg = errordlg(sprintf('Nothing was copied or the copied data is not numerical.\nThe copied data:\n\n%s', dataText));
        uiwait(dlg);
        hgclose(handles.figure1);
    end
catch
    dlg = errordlg(sprintf('Can''t read the data. Probably badly formatted (not numerical) or TOO MUCH!\nDid you copy non-sequential columns?...\nTry copying less data somehow...', dataText));
    uiwait(dlg);
    hgclose(handles.figure1);
end

% UIWAIT makes PasteDataDialog wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PasteDataDialog_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if (~isempty(handles))
    varargout{1} = handles.output;
    varargout{2} = handles.Data;
else
    varargout{1} = -1;
    varargout{2} = [];
end


% --- Executes when entered data in editable cell(s) in dataTable.
function dataTable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to dataTable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
1;


% --- Executes on button press in okButton.
function okButton_Callback(hObject, eventdata, handles)
% hObject    handle to okButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%data  = get(handles.dataTable, 'Data');
handles.Data.PastedData = get(handles.dataTable, 'Data');

columnNames = get(handles.dataTable, 'ColumnName');
if (iscell(columnNames))
    handles.Data.PastedHeaders = columnNames;
else
    handles.Data.PastedHeaders = {};
end

x = get(handles.xColPopup, 'Value');
dx = get(handles.dxColPopup, 'Value');
y = get(handles.yColPopup, 'Value');
dy = get(handles.dyColPopup, 'Value');
handles.Data.SelectedData = handles.Data.PastedData(:, [x dx y dy]);
handles.Data.WasApproved = 1;
hgclose(handles.figure1);


% --- Executes on button press in cancelButton.
function cancelButton_Callback(hObject, eventdata, handles)
% hObject    handle to cancelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Data.WasApproved = 0;
hgclose(handles.figure1);


% --- Executes on selection change in xColPopup.
function xColPopup_Callback(hObject, eventdata, handles)
% hObject    handle to xColPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns xColPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from xColPopup


% --- Executes during object creation, after setting all properties.
function xColPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xColPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in yColPopup.
function yColPopup_Callback(hObject, eventdata, handles)
% hObject    handle to yColPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns yColPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from yColPopup


% --- Executes during object creation, after setting all properties.
function yColPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yColPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in dxColPopup.
function dxColPopup_Callback(hObject, eventdata, handles)
% hObject    handle to dxColPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns dxColPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dxColPopup


% --- Executes during object creation, after setting all properties.
function dxColPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dxColPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in dyColPopup.
function dyColPopup_Callback(hObject, eventdata, handles)
% hObject    handle to dyColPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns dyColPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dyColPopup


% --- Executes during object creation, after setting all properties.
function dyColPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dyColPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
