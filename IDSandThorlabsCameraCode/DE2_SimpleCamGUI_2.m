function varargout = DE2_SimpleCamGUI_2(varargin)
% DE2_SIMPLECAMGUI MATLAB code for DE2_SimpleCamGUI.fig
%      DE2_SIMPLECAMGUI, by itself, creates a new DE2_SIMPLECAMGUI or raises the existing
%      singleton*.
%
%      H = DE2_SIMPLECAMGUI returns the handle to a new DE2_SIMPLECAMGUI or the handle to
%      the existing singleton*.
%
%      DE2_SIMPLECAMGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DE2_SIMPLECAMGUI.M with the given input arguments.
%
%      DE2_SIMPLECAMGUI('Property','Value',...) creates a new DE2_SIMPLECAMGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DE2_SimpleCamGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DE2_SimpleCamGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DE2_SimpleCamGUI

% Last Modified by GUIDE v2.5 07-Dec-2017 18:46:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DE2_SimpleCamGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @DE2_SimpleCamGUI_OutputFcn, ...
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


% --- Executes just before DE2_SimpleCamGUI is made visible.
function DE2_SimpleCamGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DE2_SimpleCamGUI (see VARARGIN)

% Choose default command line output for DE2_SimpleCamGUI
handles.output = hObject;
%Define timer properties. Timer is used to do live video
handles.tim = timer('Period', .5,                        ...
            'ExecutionMode', 'fixedSpacing',            ...
            'Timerfcn', {@timerFrame, hObject},         ...
            'Stopfcn', {@timerStop, hObject});
handles.winSz = [590 575];
handles.figure1.Position(3:4) = fliplr(handles.winSz);
handles.isDispFull.Visible = 'off';
handles.xMax = 1024;
handles.yMax = 1280;
handles.wasRun = 0;     %Flag to track when saving is on
handles.isSave = 0;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DE2_SimpleCamGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Executes when user attempts to close figure1.
% DE: Close function forces camera quit to prevent user from closing the
% program without closing the camera.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stop(handles.tim)
disp('Camera Closed')
% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Outputs from this function are returned to the command line.
function varargout = DE2_SimpleCamGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in startButton.
function startButton_Callback(hObject, eventdata, handles)
% hObject    handle to startButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Check if the time is already running. If so, do nothing. Otherwise start camera
    %This prevents user from starting two camera instances which would cause the camera to fail
if strcmp(get(handles.tim, 'Running'), 'off') && ~handles.isSave        %Only start if not already running AND not saving data
    set(handles.startButton,'string','Stop', 'backgroundcolor', 'red');
    [cam, img] = DE2_Tcam_Init(str2double(get(handles.expTime,'string')));         %Instantiate camera
    [Im, img] = DE2_Tcam_Frame(cam,img);    %Capture first frame to instantiate axes
    handles.cropIm = imagesc(Im, 'parent', handles.imAxes);           %Instantiate axes
    axes(handles.imAxes)
    axis image
    handles.fullIm = imagesc(Im, 'parent', handles.fullDisp);
    axes(handles.imAxes)
    axis image
    handles.cam = cam;                      %Add cam and img to handles for cross-function access
    handles.img = img;                      
    guidata(hObject,handles);               %Update handles prior to calling timer
    start(handles.tim);                     %Start the timer    
elseif strcmp(get(handles.tim, 'Running'), 'on') && ~handles.isSave     %Only stop if running AND not saving data
    set(handles.startButton,'string','Start', 'backgroundcolor', 'green');
    stop(handles.tim)
end
guidata(hObject, handles);

function timerFrame(~, ~, hObject)
    han = guidata(hObject);
    cam     = han.cam;
    img     = han.img;
    [Im, img] = DE2_Tcam_Frame(cam,img);
    imDat = Im;
    imDatFull = Im;
    if get(han.isDispCrop, 'value') == 1
        xC = str2double(get(han.xCrop, 'string'));
        yC = str2double(get(han.yCrop, 'string'));
        crpsz = str2double(get(han.cropSz, 'string'));
        imDat = imDat(yC-crpsz:yC+crpsz, xC-crpsz:xC+crpsz);
    end
    if (get(han.isHist, 'value') == 1) && (ishandle(han.histFig))
        set(han.histHan, 'Data', imDat);
    elseif (get(han.isHist, 'value') == 1) && ~ishandle(han.histFig)
        set(han.isHist, 'value', 0);
    end
    if get(han.isLog, 'value') ==  1
        imDat = log10(imDat+1);
        imDatFull = log10(imDatFull+1);
    end
    set(han.cropIm,'CData',imDat);
    if get(han.isDispFull, 'value') == 1
        set(han.fullIm, 'CData',imDatFull);
    end   
    if han.isSave
        fnm = [get(han.fnm, 'string') '_' num2str(get(han.tim, 'TasksExecuted')-1,'%03.0f')];   %use tasks executed to know number of executions
        if exist([fnm '.fits']) ~= 0                                                            %   -1 to index 0 to match prior saving codes
            disp('Filename already used, please change')
            set(han.saveStat, 'string', 'Error Saving')
        else
            fitswrite(imDat, fnm);
            set(han.saveStat, 'string', ['Saved ' num2str(get(han.tim, 'TasksExecuted'))])
        end                                                                 
    end    
    han.img     = img;
    set(han.minADU, 'string', num2str(min(Im(:))));
    set(han.maxADU, 'string', num2str(max(Im(:))));
    guidata(hObject, han);
    
function timerStop(~, ~, hObject)
    handles = guidata(hObject);
    cam     = handles.cam;
    DE2_Tcam_Close(cam);
    handles.isSave = 0;
    set(handles.tim,'TasksToExecute', Inf)
    set(handles.saveButton,'string','Save', 'backgroundcolor', [.47 .67 .19]);
    guidata(hObject,handles);               %Update handles prior to calling timer     
    if handles.wasRun
        handles.wasRun = 0;
        [cam, img] = DE2_Tcam_Init(str2double(get(handles.expTime,'string')));         %Instantiate camera
        handles.cam = cam;                      %Add cam and img to handles for cross-function access
        handles.img = img;                      
        guidata(hObject,handles);               %Update handles prior to calling timer
        start(handles.tim);                     %Start the timer
    end
    
    
% --- Executes on button press in ChgExp.
function ChgExp_Callback(hObject, eventdata, handles)
% hObject    handle to ChgExp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~handles.isSave      %Ignore inputs while saving
    wasRun = get(handles.tim, 'Running');
    stop(handles.tim)
    [cam, img] = DE2_Tcam_Init(str2double(get(handles.expTime,'string')));         %Instantiate camera
    handles.cam = cam;                      %Add cam and img to handles for cross-function access
    handles.img = img;                      
    guidata(hObject,handles);               %Update handles prior to calling timer
    if strcmp(wasRun, 'on')
        start(handles.tim);                     %Start the timer
    else
        DE2_Tcam_Close(cam);
    end
end

% --- Executes on button press in saveButton.
function saveButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~handles.isSave
    set(handles.saveButton,'string','Saving...', 'backgroundcolor', 'yellow');
    if strcmp(get(handles.tim, 'Running'), 'on')
        handles.wasRun = 1;
        stop(handles.tim);      %Stop current timer if present
    end
    frms = floor(str2num(get(handles.nFrames, 'string')));  %Get number of frames to record and push to integer
    set(handles.tim,'TasksToExecute', frms);     %Set execution times to nFrames
    handles.isSave  = 1;
    
    %start timer again
    [cam, img] = DE2_Tcam_Init(str2double(get(handles.expTime,'string')));         %Instantiate camera
    [Im, img] = DE2_Tcam_Frame(cam,img);    %Capture first frame to instantiate axes
    handles.cropIm = imagesc(Im, 'parent', handles.imAxes);           %Instantiate axes
    axes(handles.imAxes)
    axis image
    handles.fullIm = imagesc(Im, 'parent', handles.fullDisp);
    axes(handles.imAxes)
    axis image
    handles.cam = cam;                      %Add cam and img to handles for cross-function access
    handles.img = img;
    guidata(hObject, handles);
    start(handles.tim);                     %Start the timer  
else        %User Aborting Save
    stop(handles.tim)
end
guidata(hObject,handles);               %Update handles prior to calling timer


% --- Executes on button press in isLog.
function isLog_Callback(hObject, eventdata, handles)
% hObject    handle to isLog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of isLog

% --- Executes on button press in isDispFull.
function isDispFull_Callback(hObject, eventdata, handles)
% hObject    handle to isDispFull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject, 'Value') == 1
    handles.figure1.Position(3:4) = [960 handles.winSz(1)];
else
    handles.figure1.Position(3:4) = fliplr(handles.winSz);
end
% Hint: get(hObject,'Value') returns toggle state of isDispFull


% --- Executes on button press in isDispCrop.
function isDispCrop_Callback(hObject, eventdata, handles)
% hObject    handle to isDispCrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject, 'Value') == 1
    set(handles.isDispFull, 'visible','on')
else
    set(handles.isDispFull, 'visible','off')
    set(handles.isDispFull, 'value', 0);
    handles.figure1.Position(3:4) = fliplr(handles.winSz);
end
guidata(hObject, handles)
% Hint: get(hObject,'Value') returns toggle state of isDispCrop

% --- Executes on button press in isCrosshr.
function isCrosshr_Callback(hObject, eventdata, handles)
% hObject    handle to isCrosshr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of isCrosshr
if get(hObject, 'Value') == 1
    axes(handles.imAxes)
    handles.xCrHr = line([handles.imAxes.XLim(1) handles.imAxes.XLim(2)], [handles.imAxes.YLim(2)/2 handles.imAxes.YLim(2)/2], 'Color', 'red');
    handles.yCrHr = line([handles.imAxes.XLim(2)/2 handles.imAxes.XLim(2)/2], [handles.imAxes.YLim(1) handles.imAxes.YLim(2)], 'Color', 'red');
end
if get(hObject, 'Value') == 0
    delete(handles.xCrHr)
    delete(handles.yCrHr)
end
guidata(hObject, handles);


% --- Executes on button press in isHist.
function isHist_Callback(hObject, eventdata, handles)
% hObject    handle to isHist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of isHist
if get(hObject, 'Value') == 1
    handles.histFig = figure;
    handles.histHan = histogram(get(handles.cropIm,'CData'), 'BinLimits', [0 255], 'BinWidth', 5);
else
    if ishandle(handles.histFig)
        close(handles.histFig)
    end
end
guidata(hObject, handles)

function xCrop_Callback(hObject, eventdata, handles)
% hObject    handle to xCrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xC = str2num(get(hObject, 'string'));
crpsz = str2num(get(handles.cropSz, 'string'));
if xC - crpsz <=0
    set(handles.xCrop, 'string', num2str(crpsz + 1));
elseif xC + crpsz >= handles.xMax
    set(handles.xCrop, 'string', num2str(handles.xMax - crpsz));
end
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of xCrop as text
%        str2double(get(hObject,'String')) returns contents of xCrop as a double


% --- Executes during object creation, after setting all properties.
function xCrop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xCrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function yCrop_Callback(hObject, eventdata, handles)
% hObject    handle to yCrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yCrop as text
%        str2double(get(hObject,'String')) returns contents of yCrop as a double
yC = str2num(get(hObject, 'string'));
crpsz = str2num(get(handles.cropSz, 'string'));
if yC - crpsz <=0
    set(handles.yCrop, 'string', num2str(crpsz + 1));
elseif yC + crpsz >= handles.yMax
    set(handles.yCrop, 'string', num2str(handles.yMax - crpsz));
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function yCrop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yCrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cropSz_Callback(hObject, eventdata, handles)
% hObject    handle to cropSz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cropSz as text
%        str2double(get(hObject,'String')) returns contents of cropSz as a double
crpsz = str2num(get(hObject, 'string'));
xC = str2num(get(handles.xCrop, 'string'));
yC = str2num(get(handles.yCrop, 'string'));
if xC - crpsz <=0
    set(handles.cropSz, 'string', num2str(xC - 1));
elseif xC + crpsz > handles.xMax
    set(handles.cropSz, 'string', num2str(handles.xMax - xC));
elseif yC - crpsz <=0
    set(handles.cropSz, 'string', num2str(yC - 1));
elseif yC + crpsz > handles.yMax
    set(handles.cropSz, 'string', num2str(handles.yMax - yC));
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function cropSz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cropSz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function expTime_Callback(hObject, eventdata, handles)
% hObject    handle to expTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of expTime as text
%        str2double(get(hObject,'String')) returns contents of expTime as a double


% --- Executes during object creation, after setting all properties.
function expTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to expTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nFrames_Callback(hObject, eventdata, handles)
% hObject    handle to nFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nFrames as text
%        str2double(get(hObject,'String')) returns contents of nFrames as a double


% --- Executes during object creation, after setting all properties.
function nFrames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fnm_Callback(hObject, eventdata, handles)
% hObject    handle to fnm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fnm as text
%        str2double(get(hObject,'String')) returns contents of fnm as a double


% --- Executes during object creation, after setting all properties.
function fnm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fnm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% function varargout = DE2_SimpleCamGUI_2(varargin)
% % DE2_SIMPLECAMGUI MATLAB code for DE2_SimpleCamGUI.fig
% %      DE2_SIMPLECAMGUI, by itself, creates a new DE2_SIMPLECAMGUI or raises the existing
% %      singleton*.
% %
% %      H = DE2_SIMPLECAMGUI returns the handle to a new DE2_SIMPLECAMGUI or the handle to
% %      the existing singleton*.
% %
% %      DE2_SIMPLECAMGUI('CALLBACK',hObject,eventData,handles,...) calls the local
% %      function named CALLBACK in DE2_SIMPLECAMGUI.M with the given input arguments.
% %
% %      DE2_SIMPLECAMGUI('Property','Value',...) creates a new DE2_SIMPLECAMGUI or raises the
% %      existing singleton*.  Starting from the left, property value pairs are
% %      applied to the GUI before DE2_SimpleCamGUI_OpeningFcn gets called.  An
% %      unrecognized property name or invalid value makes property application
% %      stop.  All inputs are passed to DE2_SimpleCamGUI_OpeningFcn via varargin.
% %
% %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
% %      instance to run (singleton)".
% %
% % See also: GUIDE, GUIDATA, GUIHANDLES
% 
% % Edit the above text to modify the response to help DE2_SimpleCamGUI
% 
% % Last Modified by GUIDE v2.5 07-Dec-2017 18:18:46
% 
% % Begin initialization code - DO NOT EDIT
% gui_Singleton = 1;
% gui_State = struct('gui_Name',       mfilename, ...
%                    'gui_Singleton',  gui_Singleton, ...
%                    'gui_OpeningFcn', @DE2_SimpleCamGUI_OpeningFcn, ...
%                    'gui_OutputFcn',  @DE2_SimpleCamGUI_OutputFcn, ...
%                    'gui_LayoutFcn',  [] , ...
%                    'gui_Callback',   []);
% if nargin && ischar(varargin{1})
%     gui_State.gui_Callback = str2func(varargin{1});
% end
% 
% if nargout
%     [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
% else
%     gui_mainfcn(gui_State, varargin{:});
% end
% % End initialization code - DO NOT EDIT
% 
% 
% % --- Executes just before DE2_SimpleCamGUI is made visible.
% function DE2_SimpleCamGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% % This function has no output args, see OutputFcn.
% % hObject    handle to figure
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% % varargin   command line arguments to DE2_SimpleCamGUI (see VARARGIN)
% 
% % Choose default command line output for DE2_SimpleCamGUI
% handles.output = hObject;
% %Define timer properties. Timer is used to do live video
% handles.tim = timer('Period', .5,                        ...
%             'ExecutionMode', 'fixedSpacing',            ...
%             'Timerfcn', {@timerFrame, hObject},         ...
%             'Stopfcn', {@timerStop, hObject});
% handles.winSz = [590 575];
% handles.figure1.Position(3:4) = fliplr(handles.winSz);
% set(handles.startButton,'string','Start', 'backgroundcolor', 'green');
% handles.isDispFull.Visible = 'off';
% handles.xMax = 1024;
% handles.yMax = 1280;
% handles.isSavePress = 0;     %Flag to track when saving is on
% handles.isSave = 0;
% % Update handles structure
% guidata(hObject, handles);
% 
% % UIWAIT makes DE2_SimpleCamGUI wait for user response (see UIRESUME)
% % uiwait(handles.figure1);
% 
% % --- Executes when user attempts to close figure1.
% % DE: Close function forces camera quit to prevent user from closing the
% % program without closing the camera.
% function figure1_CloseRequestFcn(hObject, eventdata, handles)
% % hObject    handle to figure1 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% stop(handles.tim)
% disp('Camera Closed')
% % Hint: delete(hObject) closes the figure
% delete(hObject);
% 
% 
% % --- Outputs from this function are returned to the command line.
% function varargout = DE2_SimpleCamGUI_OutputFcn(hObject, eventdata, handles) 
% % varargout  cell array for returning output args (see VARARGOUT);
% % hObject    handle to figure
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Get default command line output from handles structure
% varargout{1} = handles.output;
% 
% 
% % --- Executes on button press in startButton.
% function startButton_Callback(hObject, eventdata, handles)
% % hObject    handle to startButton (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% %Check if the time is already running. If so, do nothing. Otherwise start camera
%     %This prevents user from starting two camera instances which would cause the camera to fail
% if strcmp(get(handles.tim, 'Running'), 'off')
%     set(handles.startButton,'string','Stop', 'backgroundcolor', 'red');
%     [cam, img] = DE2_Tcam_Init(str2double(get(handles.expTime,'string')));         %Instantiate camera
%     [Im, img] = DE2_Tcam_Frame(cam,img);    %Capture first frame to instantiate axes
%     handles.cropIm = imagesc(Im, 'parent', handles.imAxes);           %Instantiate axes
%     axes(handles.imAxes)
%     axis image
%     handles.fullIm = imagesc(Im, 'parent', handles.fullDisp);
%     axes(handles.imAxes)
%     axis image
%     handles.cam = cam;                      %Add cam and img to handles for cross-function access
%     handles.img = img;                      
%     guidata(hObject,handles);               %Update handles prior to calling timer
%     start(handles.tim);                     %Start the timer    
% elseif strcmp(get(handles.tim, 'Running'), 'on')
%     set(handles.startButton,'string','Start', 'backgroundcolor', 'green');
%     stop(handles.tim)
% end
% guidata(hObject, handles);
% 
% function timerFrame(~, ~, hObject)
%     han = guidata(hObject);
%     cam     = han.cam;
%     img     = han.img;
%     [Im, img] = DE2_Tcam_Frame(cam,img);
%     imDat = Im;
%     imDatFull = Im;
%     if get(han.isDispCrop, 'value') == 1
%         xC = str2double(get(han.xCrop, 'string'));
%         yC = str2double(get(han.yCrop, 'string'));
%         crpsz = str2double(get(han.cropSz, 'string'));
%         imDat = imDat(yC-crpsz:yC+crpsz, xC-crpsz:xC+crpsz);
%     end
%     if (get(han.isHist, 'value') == 1) && (ishandle(han.histFig))
%         set(han.histHan, 'Data', imDat);
%     elseif (get(han.isHist, 'value') == 1) && ~ishandle(han.histFig)
%         set(han.isHist, 'value', 0);
%     end
%     if get(han.isLog, 'value') ==  1
%         imDat = log(imDat+1);
%         imDatFull = log(imDatFull+1);
%     end
%     set(han.cropIm,'CData',imDat);
%     if get(han.isDispFull, 'value') == 1
%         set(han.fullIm, 'CData',imDatFull);
%     end        
%     han.img     = img;
%     set(han.minADU, 'string', num2str(min(Im(:))));
%     set(han.maxADU, 'string', num2str(max(Im(:))));
%     guidata(hObject, han);
%     
% function timerStop(~, ~, hObject)
%     handles = guidata(hObject);
%     cam     = handles.cam;
%     DE2_Tcam_Close(cam);
%     guidata(hObject, handles);
%     
% % --- Executes on button press in ChgExp.
% function ChgExp_Callback(hObject, eventdata, handles)
% % hObject    handle to ChgExp (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% wasRun = get(handles.tim, 'Running');
% stop(handles.tim)
% [cam, img] = DE2_Tcam_Init(str2double(get(handles.expTime,'string')));         %Instantiate camera
% handles.cam = cam;                      %Add cam and img to handles for cross-function access
% handles.img = img;                      
% guidata(hObject,handles);               %Update handles prior to calling timer
% if strcmp(wasRun, 'on')
%     start(handles.tim);                     %Start the timer
% else
%     DE2_Tcam_Close(cam);
% end
% 
% % --- Executes on button press in saveButton.
% function saveButton_Callback(hObject, eventdata, handles)
% % hObject    handle to saveButton (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% if ~handles.isSavePress
%     handles.isSavePress = 1;
% end
% 
% % --- Executes on button press in isLog.
% function isLog_Callback(hObject, eventdata, handles)
% % hObject    handle to isLog (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hint: get(hObject,'Value') returns toggle state of isLog
% 
% % --- Executes on button press in isDispFull.
% function isDispFull_Callback(hObject, eventdata, handles)
% % hObject    handle to isDispFull (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% if get(hObject, 'Value') == 1
%     handles.figure1.Position(3:4) = [960 handles.winSz(1)];
% else
%     handles.figure1.Position(3:4) = fliplr(handles.winSz);
% end
% % Hint: get(hObject,'Value') returns toggle state of isDispFull
% 
% 
% % --- Executes on button press in isDispCrop.
% function isDispCrop_Callback(hObject, eventdata, handles)
% % hObject    handle to isDispCrop (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% if get(hObject, 'Value') == 1
%     set(handles.isDispFull, 'visible','on')
% else
%     set(handles.isDispFull, 'visible','off')
%     set(handles.isDispFull, 'value', 0);
%     handles.figure1.Position(3:4) = fliplr(handles.winSz);
% end
% guidata(hObject, handles)
% % Hint: get(hObject,'Value') returns toggle state of isDispCrop
% 
% % --- Executes on button press in isCrosshr.
% function isCrosshr_Callback(hObject, eventdata, handles)
% % hObject    handle to isCrosshr (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hint: get(hObject,'Value') returns toggle state of isCrosshr
% if get(hObject, 'Value') == 1
%     axes(handles.imAxes)
%     handles.xCrHr = line([handles.imAxes.XLim(1) handles.imAxes.XLim(2)], [handles.imAxes.YLim(2)/2 handles.imAxes.YLim(2)/2], 'Color', 'red');
%     handles.yCrHr = line([handles.imAxes.XLim(2)/2 handles.imAxes.XLim(2)/2], [handles.imAxes.YLim(1) handles.imAxes.YLim(2)], 'Color', 'red');
% end
% if get(hObject, 'Value') == 0
%     delete(handles.xCrHr)
%     delete(handles.yCrHr)
% end
% guidata(hObject, handles);
% 
% 
% % --- Executes on button press in isHist.
% function isHist_Callback(hObject, eventdata, handles)
% % hObject    handle to isHist (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hint: get(hObject,'Value') returns toggle state of isHist
% if get(hObject, 'Value') == 1
%     handles.histFig = figure;
%     handles.histHan = histogram(get(handles.cropIm,'CData'), 'BinLimits', [0 255], 'BinWidth', 5);
% else
%     if ishandle(handles.histFig)
%         close(handles.histFig)
%     end
% end
% guidata(hObject, handles)
% 
% function xCrop_Callback(hObject, eventdata, handles)
% % hObject    handle to xCrop (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% xC = str2num(get(hObject, 'string'));
% crpsz = str2num(get(handles.cropSz, 'string'));
% if xC - crpsz <=0
%     set(handles.xCrop, 'string', num2str(crpsz + 1));
% elseif xC + crpsz >= handles.xMax
%     set(handles.xCrop, 'string', num2str(handles.xMax - crpsz));
% end
% guidata(hObject, handles);
% % Hints: get(hObject,'String') returns contents of xCrop as text
% %        str2double(get(hObject,'String')) returns contents of xCrop as a double
% 
% 
% % --- Executes during object creation, after setting all properties.
% function xCrop_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to xCrop (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% 
% function yCrop_Callback(hObject, eventdata, handles)
% % hObject    handle to yCrop (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of yCrop as text
% %        str2double(get(hObject,'String')) returns contents of yCrop as a double
% yC = str2num(get(hObject, 'string'));
% crpsz = str2num(get(handles.cropSz, 'string'));
% if yC - crpsz <=0
%     set(handles.yCrop, 'string', num2str(crpsz + 1));
% elseif yC + crpsz >= handles.yMax
%     set(handles.yCrop, 'string', num2str(handles.yMax - crpsz));
% end
% guidata(hObject, handles);
% 
% % --- Executes during object creation, after setting all properties.
% function yCrop_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to yCrop (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% 
% function cropSz_Callback(hObject, eventdata, handles)
% % hObject    handle to cropSz (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of cropSz as text
% %        str2double(get(hObject,'String')) returns contents of cropSz as a double
% crpsz = str2num(get(hObject, 'string'));
% xC = str2num(get(handles.xCrop, 'string'));
% yC = str2num(get(handles.yCrop, 'string'));
% if xC - crpsz <=0
%     set(handles.cropSz, 'string', num2str(xC - 1));
% elseif xC + crpsz > handles.xMax
%     set(handles.cropSz, 'string', num2str(handles.xMax - xC));
% elseif yC - crpsz <=0
%     set(handles.cropSz, 'string', num2str(yC - 1));
% elseif yC + crpsz > handles.yMax
%     set(handles.cropSz, 'string', num2str(handles.yMax - yC));
% end
% guidata(hObject, handles)
% 
% 
% % --- Executes during object creation, after setting all properties.
% function cropSz_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to cropSz (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% 
% 
% function expTime_Callback(hObject, eventdata, handles)
% % hObject    handle to expTime (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of expTime as text
% %        str2double(get(hObject,'String')) returns contents of expTime as a double
% 
% 
% % --- Executes during object creation, after setting all properties.
% function expTime_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to expTime (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% 
% function nFrames_Callback(hObject, eventdata, handles)
% % hObject    handle to nFrames (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of nFrames as text
% %        str2double(get(hObject,'String')) returns contents of nFrames as a double
% 
% 
% % --- Executes during object creation, after setting all properties.
% function nFrames_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to nFrames (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% 
% 
% 
% function fnm_Callback(hObject, eventdata, handles)
% % hObject    handle to fnm (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of fnm as text
% %        str2double(get(hObject,'String')) returns contents of fnm as a double
% 
% 
% % --- Executes during object creation, after setting all properties.
% function fnm_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to fnm (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end


% function varargout = DE2_SimpleCamGUI_2(varargin)
% % DE2_SIMPLECAMGUI MATLAB code for DE2_SimpleCamGUI.fig
% %      DE2_SIMPLECAMGUI, by itself, creates a new DE2_SIMPLECAMGUI or raises the existing
% %      singleton*.
% %
% %      H = DE2_SIMPLECAMGUI returns the handle to a new DE2_SIMPLECAMGUI or the handle to
% %      the existing singleton*.
% %
% %      DE2_SIMPLECAMGUI('CALLBACK',hObject,eventData,handles,...) calls the local
% %      function named CALLBACK in DE2_SIMPLECAMGUI.M with the given input arguments.
% %
% %      DE2_SIMPLECAMGUI('Property','Value',...) creates a new DE2_SIMPLECAMGUI or raises the
% %      existing singleton*.  Starting from the left, property value pairs are
% %      applied to the GUI before DE2_SimpleCamGUI_OpeningFcn gets called.  An
% %      unrecognized property name or invalid value makes property application
% %      stop.  All inputs are passed to DE2_SimpleCamGUI_OpeningFcn via varargin.
% %
% %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
% %      instance to run (singleton)".
% %
% % See also: GUIDE, GUIDATA, GUIHANDLES
% 
% % Edit the above text to modify the response to help DE2_SimpleCamGUI
% 
% % Last Modified by GUIDE v2.5 24-Oct-2017 13:43:53
% 
% % Begin initialization code - DO NOT EDIT
% gui_Singleton = 1;
% gui_State = struct('gui_Name',       mfilename, ...
%                    'gui_Singleton',  gui_Singleton, ...
%                    'gui_OpeningFcn', @DE2_SimpleCamGUI_OpeningFcn, ...
%                    'gui_OutputFcn',  @DE2_SimpleCamGUI_OutputFcn, ...
%                    'gui_LayoutFcn',  [] , ...
%                    'gui_Callback',   []);
% if nargin && ischar(varargin{1})
%     gui_State.gui_Callback = str2func(varargin{1});
% end
% 
% if nargout
%     [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
% else
%     gui_mainfcn(gui_State, varargin{:});
% end
% % End initialization code - DO NOT EDIT
% 
% 
% % --- Executes just before DE2_SimpleCamGUI is made visible.
% function DE2_SimpleCamGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% % This function has no output args, see OutputFcn.
% % hObject    handle to figure
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% % varargin   command line arguments to DE2_SimpleCamGUI (see VARARGIN)
% 
% % Choose default command line output for DE2_SimpleCamGUI
% handles.output = hObject;
% %Define timer properties. Timer is used to do live video
% handles.tim = timer('Period', .5,                        ...
%             'ExecutionMode', 'fixedSpacing',            ...
%             'Timerfcn', {@timerFrame, hObject},         ...
%             'Stopfcn', {@timerStop, hObject});
% 
% % Update handles structure
% guidata(hObject, handles);
% 
% % UIWAIT makes DE2_SimpleCamGUI wait for user response (see UIRESUME)
% % uiwait(handles.figure1);
% 
% 
% % --- Outputs from this function are returned to the command line.
% function varargout = DE2_SimpleCamGUI_OutputFcn(hObject, eventdata, handles) 
% % varargout  cell array for returning output args (see VARARGOUT);
% % hObject    handle to figure
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Get default command line output from handles structure
% varargout{1} = handles.output;
% 
% 
% % --- Executes on button press in startButton.
% function startButton_Callback(hObject, eventdata, handles)
% % hObject    handle to startButton (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% %Check if the time is already running. If so, do nothing. Otherwise start camera
%     %This prevents user from starting two camera instances which would cause the camera to fail
% if strcmp(get(handles.tim, 'Running'), 'off')
%     set(handles.startButton,'string','Stop', 'backgroundcolor', 'red');
%     [cam, img] = DE2_Tcam_Init(.1);         %Instantiate camera
%     [Im, img] = DE2_Tcam_Frame(cam,img);    %Capture first frame to instantiate axes
%     handles.imAxis = imagesc(Im);           %Instantiate axes
%     handles.cam = cam;                      %Add cam and img to handles for cross-function access
%     handles.img = img;                      
%     guidata(hObject,handles);               %Update handles prior to calling timer
%     start(handles.tim);                     %Start the timer    
% elseif strcmp(get(handles.tim, 'Running'), 'on')
%     set(handles.startButton,'string','Start', 'backgroundcolor', 'green');
%     stop(handles.tim)
% end
% guidata(hObject, handles);
% 
% 
% % --- Executes on button press in StopButton.
% function StopButton_Callback(hObject, eventdata, handles)
% % hObject    handle to StopButton (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% stop(handles.tim);
% %handles.camflg = 0;
% guidata(hObject, handles);
% 
% function timerFrame(~, ~, hObject)
%     handles = guidata(hObject);
%     cam     = handles.cam;
%     img     = handles.img;
%     [Im, img] = DE2_Tcam_Frame(cam,img);
%     set(handles.imAxis,'CData',Im);
%     handles.img     = img;
%     disp('t');
%     guidata(hObject, handles);
%     
% function timerStop(~, ~, hObject)
%     handles = guidata(hObject);
%     cam     = handles.cam;
%     DE2_Tcam_Close(cam);
%     guidata(hObject, handles);
