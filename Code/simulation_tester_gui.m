function varargout = simulation_tester_gui(varargin)
% SIMULATION_TESTER_GUI MATLAB code for simulation_tester_gui.fig
%      SIMULATION_TESTER_GUI, by itself, creates a new SIMULATION_TESTER_GUI or raises the existing
%      singleton*.
%
%      H = SIMULATION_TESTER_GUI returns the handle to a new SIMULATION_TESTER_GUI or the handle to
%      the existing singleton*.
%
%      SIMULATION_TESTER_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SIMULATION_TESTER_GUI.M with the given input arguments.
%
%      SIMULATION_TESTER_GUI('Property','Value',...) creates a new SIMULATION_TESTER_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before simulation_tester_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to simulation_tester_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help simulation_tester_gui

% Last Modified by GUIDE v2.5 20-Apr-2015 20:28:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @simulation_tester_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @simulation_tester_gui_OutputFcn, ...
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

% --- Executes just before simulation_tester_gui is made visible.
function simulation_tester_gui_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to simulation_tester_gui (see VARARGIN)

% Choose default command line output for simulation_tester_gui
handles.output = hObject;
handles.data_array = zeros(7,2);
handles.degr_min = 150;
handles.freq = 2;
handles.amp = 10;
handles.num = 200;
handles.length = 50;
handles.max_freq = 5;
handles.line = 0;
handles.min_range = 0;
handles.max_range = 20;
handles.n_range = 100;
handles.vessel_type_num = 0;
handles.xl = 0;
handles.yl = 0;
% Values specific to experimental 
handles.voxel_dimxy = 0.375*10^-3;
handles.voxel_dimz = 1*10^-3;
handles.spacing = 0.5*10^-3;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes simulation_tester_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = simulation_tester_gui_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, ~, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
cla(handles.sim_fig)
clearvars x y z contents selected u v N R

contents = cellstr(get(hObject,'String'));
selected = contents{get(hObject,'Value')};

v = linspace(0,1,handles.num);

A = handles.amp;

handles.vessel_type_num = str2double(selected(1));

% Vessel types
if str2double(selected(1)) == 1 % Straight
    x = A*ones(1,handles.num);
    y = A*ones(1,handles.num);
    z = v*handles.length;
elseif str2double(selected(1)) == 2 % Sinusoidal
    x=A*sin(2*pi*v*handles.freq);
    y=ones(1,handles.num);
    z=v*handles.length;
elseif str2double(selected(1)) == 3 % Coiled
    x=A*cos(2*pi*v*handles.freq);
    y=A*sin(2*pi*v*handles.freq);
    z=v*handles.length;
elseif str2double(selected(1)) == 4 % Highly Coiled
    x=A*cos(2*pi*v*handles.freq);
    y=A*sin(2*pi*v*handles.freq)+A*cos(5*pi*v*handles.freq);
    z=v*handles.length;
end

plot3(x,y,z,'Parent',handles.sim_fig,'Marker','.','LineStyle',':'); %plot centreline
grid on

sim_cent = zeros(size(x,2),3);
sim_cent(:,1)=x;
sim_cent(:,2)=y;
sim_cent(:,3)=z;

handles.data_array(:,1) = centreline_tort_beta(sim_cent,[handles.voxel_dimxy, ...
    handles.voxel_dimxy,handles.voxel_dimz],handles.spacing);
assignin('base','d',handles.data_array(5))
assignin('base','centre',sim_cent)
set(handles.tort_table_edit,'Data',handles.data_array)

guidata(hObject, handles);
    

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject,~, ~)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function length_text_edit_Callback(~, ~, ~)
% hObject    handle to length_text_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of length_text_edit as text
%        str2double(get(hObject,'String')) returns contents of length_text_edit as a double

% --- Executes during object creation, after setting all properties.
function length_text_edit_CreateFcn(hObject, ~, ~)
% hObject    handle to length_text_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function freq_edit_Callback(hObject, ~, handles)
% hObject    handle to freq_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of freq_edit as text
%        str2double(get(hObject,'String')) returns contents of freq_edit as a double
handles.freq = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function freq_edit_CreateFcn(hObject, ~, ~)
% hObject    handle to freq_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function amp_edit_Callback(hObject, ~, handles)
% hObject    handle to amp_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of amp_edit as text
%        str2double(get(hObject,'String')) returns contents of amp_edit as a double
handles.amp = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function amp_edit_CreateFcn(hObject, ~, ~)
% hObject    handle to amp_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function num_edit_Callback(hObject, ~, handles)
% hObject    handle to num_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_edit as text
%        str2double(get(hObject,'String')) returns contents of num_edit as a double
handles.num=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function num_edit_CreateFcn(hObject, ~, ~)
% hObject    handle to num_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function deg_min_edit_Callback(hObject, ~, handles)
% hObject    handle to deg_min_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of deg_min_edit as text
%        str2double(get(hObject,'String')) returns contents of deg_min_edit as a double
handles.degr_min = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function deg_min_edit_CreateFcn(hObject, ~, ~)
% hObject    handle to deg_min_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function length_edit_Callback(hObject, ~, handles)
% hObject    handle to length_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of length_edit as text
%        str2double(get(hObject,'String')) returns contents of length_edit as a double
handles.length=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function length_edit_CreateFcn(hObject, ~, ~)
% hObject    handle to length_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in plot_freq.
function plot_freq_Callback(~, ~,~)
% hObject    handle to plot_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function max_freq_Callback(hObject, ~, handles)
% hObject    handle to max_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_freq as text
%        str2double(get(hObject,'String')) returns contents of max_freq as a double
handles.max_freq = str2double(get(hObject,'String'));

guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function max_freq_CreateFcn(hObject, ~, ~)
% hObject    handle to max_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function min_edit_Callback(hObject, ~, handles)
% hObject    handle to min_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_edit as text
%        str2double(get(hObject,'String')) returns contents of min_edit as a double
handles.min_range = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function min_edit_CreateFcn(hObject, ~, ~)
% hObject    handle to min_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function max_edit_Callback(hObject, ~, handles)
% hObject    handle to max_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_edit as text
%        str2double(get(hObject,'String')) returns contents of max_edit as a double
handles.max_range = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function max_edit_CreateFcn(hObject, ~, ~)
% hObject    handle to max_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in variable_range_select.
function variable_range_select_Callback(hObject, ~, handles)
% hObject    handle to variable_range_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns variable_range_select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from variable_range_select
cla(handles.range_fig)

contents = cellstr(get(hObject,'String'));
selected = contents{get(hObject,'Value')};

% ------------------------------
% Parameter Modulation 
% ------------------------------
v = linspace(0,1,handles.num);
line = zeros(handles.n_range,handles.num,3);
range = linspace(handles.min_range,handles.max_range,handles.n_range);
A = handles.amp;

% ---------
% Frequency
% ---------
if str2double(selected(1)) == 1 % Frequency
    % testing range from min->max with n intervals
    if handles.vessel_type_num == 1 % Straight (no change)
        x = A*ones(1,handles.num);
        y = A*ones(1,handles.num);
        z = v*handles.length;
        for n=1:handles.n_range;        
            line(n,:,1) = x;
            line(n,:,2) = y;
            line(n,:,3) = z;
        end
    elseif handles.vessel_type_num == 2 % Sine (x->2D)
        x=zeros(handles.n_range,handles.num);
        x(:,:)=A*sin(2*pi*range'*v);
        y=ones(1,handles.num);
        z=v*handles.length;
        for n=1:handles.n_range;
            line(n,:,1)=x(n,:);
            line(n,:,2)=y;
            line(n,:,3)=z;
        end
    elseif handles.vessel_type_num == 3 % Coil (x&y->2D)
        x=zeros(handles.n_range,handles.num);
        y=zeros(handles.n_range,handles.num);
        x(:,:)=A*cos(2*pi*range'*v);
        y(:,:)=A*sin(2*pi*range'*v);
        z=v*handles.length;
        for n=1:handles.n_range;
            line(n,:,1)=x(n,:);
            line(n,:,2)=y(n,:);
            line(n,:,3)=z;
        end
    elseif handles.vessel_type_num == 4 % V-Coiled (x&y->2D)
        x=zeros(handles.n_range,handles.num);
        y=zeros(handles.n_range,handles.num);
        x(:,:)=A*cos(2*pi*range'*v);
        y(:,:)=A*sin(2*pi*range'*v)+A*cos(5*pi*range'*v);
        z=v*handles.length;
        for n=1:handles.n_range;
            line(n,:,1)=x(n,:);
            line(n,:,2)=y(n,:);
            line(n,:,3)=z;
        end
    end
    
%----------
%Amplitude
%----------
elseif str2double(selected(1)) == 2
    if handles.vessel_type_num == 1 % Straight (x&y->2D)
        x=zeros(handles.n_range,handles.num);
        y=zeros(handles.n_range,handles.num);
        x(:,:) = range'*ones(1,handles.num);
        y(:,:) = range'*ones(1,handles.num);
        z = v*handles.length;
        for n=1:handles.n_range;
            line(n,:,1)=x(n,:);
            line(n,:,2)=y(n,:);
            line(n,:,3)=z;
        end
    elseif handles.vessel_type_num == 2 % Sine (x&y->2D)
        x=zeros(handles.n_range,handles.num);
        y=zeros(handles.n_range,handles.num);
        x(:,:)=range'*sin(2*pi*v*handles.freq);
        y(:,:)=range'*ones(1,handles.num);
        z=v*handles.length;
        for n=1:handles.n_range;
            line(n,:,1)=x(n,:);
            line(n,:,2)=y(n,:);
            line(n,:,3)=z;
        end
    elseif handles.vessel_type_num == 3 % Coil (x&y->2D)
        x=zeros(handles.n_range,handles.num);
        y=zeros(handles.n_range,handles.num);
        x(:,:)=range'*cos(2*pi*v*handles.freq);
        y(:,:)=range'*sin(2*pi*v*handles.freq);
        z=v*handles.length;
        for n=1:handles.n_range;
            line(n,:,1)=x(n,:);
            line(n,:,2)=y(n,:);
            line(n,:,3)=z;
        end
    elseif handles.vessel_type_num == 4 % Highly Coiled (x&y->2D)
        x=zeros(handles.n_range,handles.num);
        y=zeros(handles.n_range,handles.num);
        x(:,:)=range'*cos(2*pi*v*handles.freq);
        y(:,:)=range'*sin(2*pi*v*handles.freq)+range'*cos(5*pi*v*handles.freq);
        z=v*handles.length;
        for n=1:handles.n_range;
            line(n,:,1)=x(n,:);
            line(n,:,2)=y(n,:);
            line(n,:,3)=z;
        end
    end
    
%------------------
% Number of Points
%------------------
% number of points parameter here gives us number of steps 
% e.g. 1 point:200 points in 100 steps
elseif str2double(selected(1)) == 3
    v = zeros(handles.n_range,handles.max_range);
    for n = 1:handles.max_range
        v(n,1:n)= linspace(0,1,n);
    end

    if handles.vessel_type_num == 1 % Straight
        x = A*v;
        y = A*v;
        z = v*handles.length;
        for n=1:handles.max_range;
            line(n,1:n,1)=x(n,1:n);
            line(n,1:n,2)=y(n,1:n);
            line(n,1:n,3)=z(n,1:n);
        end
    elseif handles.vessel_type_num == 2 % Sine
        x=A*sin(2*pi*v*handles.freq);
        y=A*v;
        z=v*handles.length;
        for n=1:handles.max_range;
            line(n,1:n,1)=x(n,1:n);
            line(n,1:n,2)=y(n,1:n);
            line(n,1:n,3)=z(n,1:n);
        end
    elseif handles.vessel_type_num == 3 % Coil
        x=A*cos(2*pi*v*handles.freq);
        y=A*sin(2*pi*v*handles.freq);
        z=v*handles.length;
        for n=1:handles.max_range;
            line(n,1:n,1)=x(n,1:n);
            line(n,1:n,2)=y(n,1:n);
            line(n,1:n,3)=z(n,1:n);
        end
    elseif handles.vessel_type_num == 4 % Highly Coiled
         x=A*cos(2*pi*v*handles.freq);
         y=A*sin(2*pi*v*handles.freq)+A*cos(5*pi*v*handles.freq);
         z=v*handles.length;
         for n=1:handles.max_range;
            line(n,1:n,1)=x(n,1:n);
            line(n,1:n,2)=y(n,1:n);
            line(n,1:n,3)=z(n,1:n);
         end
    end


% --------
% Length
% --------
elseif str2double(selected(1)) == 4 
    if handles.vessel_type_num == 1 % Straight (z->2D)
       z=zeros(handles.n_range,handles.num);
       x = A*ones(1,handles.num);
       y = A*ones(1,handles.num);
       z(:,:) = range'*v;
       for n=1:handles.n_range;
           line(n,:,1)=x;
           line(n,:,2)=y;
           line(n,:,3)=z(n,:);
       end
    elseif handles.vessel_type_num == 2 % Sine (z->2D)
       z=zeros(handles.n_range,handles.num);
       x=A*sin(2*pi*v*handles.freq);
       y=ones(1,handles.num);
       z(:,:)=range'*v;
       for n=1:handles.n_range;
           line(n,:,1)=x;
           line(n,:,2)=y;
           line(n,:,3)=z(n,:);
       end
    elseif handles.vessel_type_num == 3 % Coil (z->2D)
        z=zeros(handles.n_range,handles.num);
        x=A*cos(2*pi*v*handles.freq);
        y=A*sin(2*pi*v*handles.freq);
        z(:,:)=range'*v;
        for n=1:handles.n_range;
            line(n,:,1)=x;
            line(n,:,2)=y;
            line(n,:,3)=z(n,:);
        end
    elseif handles.vessel_type_num == 4 % Highly Coiled (z->2D)
        z=zeros(handles.n_range,handles.num);        
        x=A*cos(2*pi*v*handles.freq);
        y=A*sin(2*pi*v*handles.freq)+A*cos(5*pi*v*handles.freq);
        z(:,:)=range'*v;
        for n=1:handles.n_range;
            line(n,:,1)=x;
            line(n,:,2)=y;
            line(n,:,3)=z(n,:);
        end
    end
end


assignin('base','line',line)

% Calculate metrics for ranged centreline
data_values = zeros(handles.n_range,7);
if str2double(selected(1)) == 3 % If number of points paramater test
    for nop = 2:handles.n_range
        data_values(nop,:) = centreline_tort_beta(squeeze(line(1:nop,1:nop,:)),[handles.voxel_dimxy, ...
            handles.voxel_dimxy,handles.voxel_dimz],handles.spacing);
    end
    data_values(isnan(data_values)) = 0;
    data_values(isinf(data_values)) = 0;
else 
    for n = 1:handles.n_range
        data_values(n,:) = centreline_tort_beta(squeeze(line(n,:,:)),[handles.voxel_dimxy, ...
            handles.voxel_dimxy,handles.voxel_dimz],handles.spacing);
    end
end

% Axes of ranged modulation
axes_dat = linspace(handles.min_range,handles.max_range,handles.n_range);

assignin('base','data',data_values)



% Calculate errors in values
[fitobject1,gof1]=fit(axes_dat',squeeze(data_values(:,1)),'smoothingspline'); %#ok<ASGLU>
[fitobject2,gof2]=fit(axes_dat',squeeze(data_values(:,2)),'smoothingspline'); %#ok<ASGLU>
[fitobject3,gof3]=fit(axes_dat',squeeze(data_values(:,3)),'smoothingspline'); %#ok<ASGLU>
[fitobject4,gof4]=fit(axes_dat',squeeze(data_values(:,4)),'smoothingspline'); %#ok<ASGLU>
[fitobject6,gof6]=fit(axes_dat',squeeze(data_values(:,6)),'smoothingspline'); %#ok<ASGLU>
[fitobject7,gof7]=fit(axes_dat',squeeze(data_values(:,7)),'smoothingspline'); %#ok<ASGLU>

% Extract standard error
handles.data_array(1,2)=gof1.rmse;
handles.data_array(2,2)=gof2.rmse;
handles.data_array(3,2)=gof3.rmse;
handles.data_array(4,2)=gof4.rmse;
handles.data_array(6,2)=gof6.rmse;
handles.data_array(7,2)=gof7.rmse;

% Update data table on gui
set(handles.tort_table_edit,'Data',handles.data_array)

% Plot data
set(handles.range_fig, 'YLim',[0, round(max(max(data_values)))])
hold(handles.range_fig, 'on' )
plot(axes_dat,data_values(:,1),'Parent',handles.range_fig,'Color','r','Marker','+','LineStyle','none')
plot(axes_dat,data_values(:,2),'Parent',handles.range_fig,'Color','b','Marker','*','LineStyle','none')
plot(axes_dat,data_values(:,3),'Parent',handles.range_fig,'Color','g','Marker','o','LineStyle','none')
plot(axes_dat,data_values(:,4),'Parent',handles.range_fig,'Color','k','Marker','s','LineStyle','none')
plot(axes_dat,data_values(:,6),'Parent',handles.range_fig,'Color','y','Marker','d','LineStyle','none')
plot(axes_dat,data_values(:,7),'Parent',handles.range_fig,'Color','m','Marker','x','LineStyle','none')
hold(handles.range_fig,'off')
legend(handles.range_fig,'DM','SOAM','ICMn','ICMb','StdAvdK','NormK');


handles.xl = xlim;
handles.yl = ylim;

set(handles.y_slid,'Max',handles.yl(2))
set(handles.x_slid,'Max',handles.xl(2))
set(handles.y_slid,'Value',handles.yl(2))
set(handles.x_slid,'Value',handles.xl(2))

guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function variable_range_select_CreateFcn(hObject,~,~)
% hObject    handle to variable_range_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function n_range_Callback(hObject, ~, handles)
% hObject    handle to n_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n_range as text
%        str2double(get(hObject,'String')) returns contents of n_range as a double
handles.n_range = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function n_range_CreateFcn(hObject, ~, ~)
% hObject    handle to n_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function y_slid_Callback(hObject, ~, handles)
% hObject    handle to y_slid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

set(hObject,'Max',100)
set(hObject,'Min',0)
y_max = get(hObject,'Value');

ylim(handles.range_fig, [0,y_max]);

guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function y_slid_CreateFcn(hObject, ~, ~)
% hObject    handle to y_slid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function x_slid_Callback(hObject, ~, handles)
% hObject    handle to x_slid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

x_max = get(hObject,'Value');
xlim(handles.range_fig, [0,x_max]);

guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function x_slid_CreateFcn(hObject, ~, ~)
% hObject    handle to x_slid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
