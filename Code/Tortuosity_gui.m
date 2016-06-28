function varargout = Tortuosity_gui(varargin)
%TORTUOSITY_GUI M-file for Tortuosity_gui.fig
%      TORTUOSITY_GUI, by itself, creates a new TORTUOSITY_GUI or raises the existing
%      singleton*.
%
%      H = TORTUOSITY_GUI returns the handle to a new TORTUOSITY_GUI or the handle to
%      the existing singleton*.
%
%      TORTUOSITY_GUI('Property','Value',...) creates a new TORTUOSITY_GUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to Tortuosity_gui_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      TORTUOSITY_GUI('CALLBACK') and TORTUOSITY_GUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in TORTUOSITY_GUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Tortuosity_gui

% Last Modified by GUIDE v2.5 17-Apr-2015 19:48:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Tortuosity_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @Tortuosity_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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

% Define structure to hold variables


% --- Executes just before Tortuosity_gui is made visible.
function Tortuosity_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for Tortuosity_gui

handles.data = zeros(512,512,512);
handles.data_thresh_smooth = zeros(512,512,512);
handles.data_table = zeros(6,3);
handles.XL = [0 512];    % X Zoom Limit
handles.YL = [0 512];    % Y Zoom Limit
handles.output = hObject;
handles.voxel_dim = zeros(1,3);
handles.spacing = 0;
handles.fixed = 0;
handles.x_clean = 0;
handles.y_clean = 0;
handles.z_clean = 0;
handles.thresh = 350;      % Threshold Intensity
handles.slice  = 10;     % Image Z - Slice
handles.X_size = 10;
handles.Y_size = 10;
handles.Z_size = 10;
handles.max_I = 1500;
handles.X_min = 1;
handles.Y_min = 1;
handles.Z_min = 1;
handles.X_max = 1;
handles.Y_max = 1;
handles.Z_max = 1;
handles.length = 0;
handles.dist_metric = 0; 
handles.vert_x = 0;
handles.av_proj_ang =0;
handles.largestAllowableDistance = 8;
handles.numAllowableNeighbors = 11;
handles.degr_min = 165;             % Minimum inflection inversion angle
handles.skel = 0;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Tortuosity_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Tortuosity_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in centreline_button.
function centreline_button_Callback(hObject, eventdata, handles)
% hObject    handle to centreline_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

XL=round(handles.XL);
YL=round(handles.YL);

% Skeletization of vessel to extract centreline (Skeleton3d)
start = handles.data_thresh_smooth(YL(1):YL(2),XL(1):XL(2),handles.Z_min:handles.Z_max);

% Make the image binary
start(start > 0 ) = 1;

% Remove empty voxels within coordinate limits
% -------------  UNKNOWN  ------------------

% Find centreline through thinning algorithm
handles.skel  = Skeleton3D(start);

% Extract centreline coordinate locations
w=size(handles.skel,1);
l=size(handles.skel,2);
h=size(handles.skel,3);
[x,y,z]=ind2sub([w,l,h],find(handles.skel(:)));
handles.x_clean = x;
handles.y_clean = y;
handles.z_clean = z;

% Remove noise by nearest neighbour distance comparison
fixed_xyz = [x,y,z]; % Initialize "cleaned" output array.
assignin('base','origin','fixed_xyz')
rowsToDelete = [];
for k = 1 : length(x);
    thisX = x(k);
    thisY = y(k);
    thisZ = z(k);
    distances = sqrt((thisX-x).^2 + (thisY-y).^2 + (thisZ - z).^2);
    [sortedDistances, ~] = sort(distances, 'ascend');
    if sortedDistances(handles.numAllowableNeighbors) > handles.largestAllowableDistance
        % If outlier - record row of noise
		rowsToDelete  = [rowsToDelete , k];
	else
        % It's near the curve.
    end
end

% Remove the rows of noise
fixed_xyz(rowsToDelete, :) = [];
assignin('base','origin2','fixed_xyz')

axes(handles.clean_midline_plot);
% Plot cleaned centreline_button
plot3(fixed_xyz(:,2),fixed_xyz(:,1),fixed_xyz(:,3),...
    'square','Markersize',2,'MarkerFaceColor','b','Color','b',...
    'Parent',handles.clean_midline_plot);
view(handles.clean_midline_plot,[-39,6])
grid on
% beyond this brush outliers with push button tool
handles.fixed = fixed_xyz;

% Calculate parameters of centreline
handles.data_table(:,1)=centreline_tort(fixed_xyz,handles.voxel_dim);

% Load data into table
set(handles.tort_table,'Data',handles.data_table)

guidata(hObject, handles);
   


function thresh_set_Callback(hObject, eventdata, handles) %#ok<I%#ok<MSNU> NUSL,DEFNU>
% hObject    handle to thresh_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresh_set as text
%        str2double(get(hObject,'String')) returns contents of thresh_set as a double

% Get threshold from edit box
handles.thresh = str2double(get(hObject,'String'));

% Create temp array data_thresh to store thresholded data
data_thresh = handles.data;
data_thresh(data_thresh<handles.thresh) = 0;

% Plot thresholded image
imshow(squeeze(data_thresh(:,:,handles.slice)));
axes(handles.initial_plot)
colormap(jet)
colorbar
caxis([0 handles.max_I])
xlim([handles.XL])
ylim([handles.YL])
drawnow;

clearvars  data_thresh

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function thresh_set_CreateFcn(hObject, eventdata, handles) %#ok<INUSD,DE%#ok<MSNU> FNU>
% hObject    handle to thresh_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function slice_slider_Callback(hObject, eventdata, handles) %#ok<*INUSL,DEF%#ok<MSNU> NU>
% hObject    handle to slice_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Set and get slice values
b = get(hObject,'Value');
set(hObject,'Value',b)
handles.slice = round(b);

% Set limits of axes and thresholds
handles.XL = xlim(handles.initial_plot);
handles.YL = ylim(handles.initial_plot);
data_thresh = handles.data;
data_thresh(data_thresh<handles.thresh) = 0;
data_thresh = squeeze(data_thresh(:,:,handles.slice));

% Plot slices
imshow(data_thresh,'Parent',handles.initial_plot);
axes(handles.initial_plot)
colormap(jet)
colorbar
caxis(handles.initial_plot, [0 handles.max_I])
xlim([handles.XL])
ylim([handles.YL])

clearvars data_thresh

guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function slice_slider_CreateFcn(hObject, eventdata, handles) %#ok<INUSD,*DEFNU>
% hObject    handle to slice_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in open_files_pshbut.
function open_files_pshbut_Callback(hObject, eventdata, handles)
% hObject    handle to open_files_pshbut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Set current folder
cd('C:\Users\Will\Desktop') 
[filename,pathname] = uigetfile('*.*','MultiSelect','on');

ext = regexp(filename, '(?<=\.)\d+$', 'match', 'once');
[~, order] = sort(str2double(ext));
filename = filename(order);

set(handles.patient_num,'String', pathname)

% Find data array size
s = dicomread([pathname,char(filename(1))]);
a = dicominfo([pathname,char(filename(1))]);
handles.X_size = size(s,1);
handles.Y_size = size(s,2);
handles.Z_size = numel(filename);

% Voxel dimensions
handles.voxel_dim(1) = a.PixelSpacing(1);
handles.voxel_dim(2) = a.PixelSpacing(2);
handles.voxel_dim(3) = a.SliceThickness;
handles.spacing      = a.SpacingBetweenSlices;
% Load data into array 'X'
handles.data = zeros(handles.X_size,handles.Y_size,handles.Z_size);
for n = 1:handles.Z_size
    handles.data(:,:,n) = dicomread([pathname,char(filename(n))]);
end

% Plot data
imshow(squeeze(handles.data(:,:,handles.slice)),'Parent',handles.initial_plot);
axes(handles.initial_plot)
colormap(jet)
colorbar
caxis(handles.initial_plot, [0 handles.max_I])

% set slider slice values
set(handles.slice_slider,'Min',1);
set(handles.slice_slider,'Max',handles.Z_size);
set(handles.slice_slider,'SliderStep',[(1/handles.Z_size) , (10/handles.Z_size) ]);

clearvars filename pathname


guidata(hObject,handles)

% --- Executes on button press in select_roi_pshbut.
function select_roi_pshbut_Callback(hObject, eventdata, handles)
% hObject    handle to select_roi_pshbut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get axes limits
handles.XL = xlim(handles.initial_plot);
handles.YL = ylim(handles.initial_plot);

% Isolate data
handles.data_thresh_smooth = handles.data;
handles.data_thresh_smooth(handles.data_thresh_smooth<handles.thresh) = 0;
imshow(squeeze(handles.data_thresh_smooth(:,:,handles.slice)),...
    'Parent',handles.initial_plot)

% Zoom to predefined limits & Set Color
xlim(handles.XL)
ylim(handles.YL)
axes(handles.initial_plot)
colormap(jet)
colorbar
caxis(handles.initial_plot, [0 handles.max_I])

% Create R.O.I. mask
mask = roipoly;

% Load data with ROI
roi = handles.data_thresh_smooth; 

% Use binary mask to isolate ROI data [0/1]
for x=1:size(mask,1)
    for y = 1:size(mask,2)
        if mask(x,y) == 0
            roi(x,y,:) = 0;
        end
    end
end

% Max size of roi
handles.X_max = size(handles.data_thresh_smooth,1);
handles.Y_max = size(handles.data_thresh_smooth,2);
handles.Z_max = size(handles.data_thresh_smooth,3);
handles.Z_min=1;

handles.data_thresh_smooth =roi;
clearvars data_thresh mask

% Plot ROI data points
 cla(handles.roi_axes)

% Plot isosurface
[xs,ys,zs,D] = subvolume(roi,...
            [handles.X_min,handles.X_max,...
            handles.Y_min,handles.Y_max,...
            handles.Z_min,handles.Z_max]); %Isolate subvolume

xlabel(handles.roi_axes,'x')
ylabel(handles.roi_axes,'y')
zlabel(handles.roi_axes,'z')
axes(handles.roi_axes);
patch(isosurface(xs,ys,zs,D),'FaceColor','red','EdgeColor','none','Parent',handles.roi_axes)
view(handles.roi_axes,[-45,16])
set(handles.roi_axes,'YTickLabel',[],'XTickLabel',[])
grid on
camlight
lighting gouraud

clearvars x y z D roi

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function roi_axes_CreateFcn(hObject, eventdata, handles) %#ok<INUSD>
% hObject    handle to roi_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate roi_axes



function max_intensity_Callback(hObject, eventdata, handles)
% hObject    handle to max_intensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_intensity as text
%        str2double(get(hObject,'String')) returns contents of max_intensity as a double

% Get and set maximum image intensity
handles.max_I = str2double(get(hObject, 'String'));
caxis(handles.initial_plot,[0 handles.max_I])
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function max_intensity_CreateFcn(hObject, eventdata, handles) %#ok<INUSD>
% hObject    handle to max_intensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function Zmaxroi_Callback(hObject, eventdata, handles)
% hObject    handle to Zmaxroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% get and set maximum Z value in the ROI plot
handles.Z_max = round(handles.Z_size - get(hObject,'Value'));
zlim(handles.roi_axes, [handles.Z_min,handles.Z_max]);

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function Zmaxroi_CreateFcn(hObject, eventdata, handles) %#ok<INUSD>
% hObject    handle to Zmaxroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function Zminroi_Callback(hObject, eventdata, handles)
% hObject    handle to Zminroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% get and set Z min in ROI plot
handles.Z_min = round(1 + get(hObject,'Value'));
zlim(handles.roi_axes, [handles.Z_min, handles.Z_max])

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Zminroi_CreateFcn(hObject, eventdata, handles) %#ok<INUSD>
% hObject    handle to Zminroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function nearest_neigh_slid_Callback(hObject, eventdata, handles)
% hObject    handle to nearest_neigh_slid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Get number of nearest neighbour limit
new_neigh = round(get(hObject,'Value'));
set(hObject,'Value',new_neigh)

% Remove noise voxels by distance comparisons
x = handles.x_clean;
y = handles.y_clean;
z = handles.z_clean;
% Initialize "cleaned" output array.
fixed_xyz = [x, y, z];
% Initialize location of noise array 
rowsToDelete = [];
for k = 1 : length(x);
    thisX = x(k);
    thisY = y(k);
    thisZ = z(k);
    distances = sqrt((thisX-x).^2 + (thisY-y).^2 + (thisZ - z).^2);
    [sortedDistances, ~] = sort(distances, 'ascend');
    if sortedDistances(new_neigh) > handles.largestAllowableDistance
        % If outlier - record row of noise
		rowsToDelete  = [rowsToDelete , k]; %#ok<*AGROW>
	else
        % It's near the curve.
    end
end
% Remove the rows.
fixed_xyz(rowsToDelete, :) = [];

% Plot cleaned centreline_button
plot3(fixed_xyz(:,2),fixed_xyz(:,1),fixed_xyz(:,3),...
    'square','Markersize',2,'MarkerFaceColor','b','Color','b',...
    'Parent',handles.clean_midline_plot);
view(handles.clean_midline_plot,[44,6])
zlim(handles.clean_midline_plot,[0 max(z)])

guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function nearest_neigh_slid_CreateFcn(hObject, eventdata, handles) %#ok<INUSD>
% hObject    handle to nearest_neigh_slid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on slider movement.
function min_dist_slid_Callback(hObject, eventdata, handles)
% hObject    handle to min_dist_slid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% get and set minimum distance value
new_d = round(get(hObject,'Value'));
set(hObject,'Value',new_d);

% remove noise by distance comparisons
x = handles.x_clean;
y = handles.y_clean;
z = handles.z_clean;
fixed_xyz = [x, y, z]; % Initialize our "cleaned" output array.
rowsToDelete = [];
for k = 1 : length(x);
    thisX = x(k);
    thisY = y(k);
    thisZ = z(k);
    distances = sqrt((thisX-x).^2 + (thisY-y).^2 + (thisZ - z).^2);
    [sortedDistances, ~] = sort(distances, 'ascend');
    if sortedDistances(handles.numAllowableNeighbors) > new_d
        % If outlier - record row of noise
		rowsToDelete  = [rowsToDelete , k];
	else
        % It's near the curve.
    end
end
% Remove the rows.
fixed_xyz(rowsToDelete, :) = [];

% Plot cleaned centreline_button
plot3(fixed_xyz(:,2),fixed_xyz(:,1),fixed_xyz(:,3),...
    'square','Markersize',2,'MarkerFaceColor','b','Color','b',...
    'Parent',handles.clean_midline_plot);
view(handles.clean_midline_plot,[44,6])
zlim(handles.clean_midline_plot,[0 max(z)])

guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function min_dist_slid_CreateFcn(hObject, eventdata, handles) %#ok<INUSD>
% hObject    handle to min_dist_slid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Brush outlier tool
function uitoggletool8_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

brush('on');
   
guidata(hObject, handles)


% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Extract data from figure

% handles to low-level graphics object
dataObjs = get(handles.clean_midline_plot, 'Children'); 
% low-level grahics object to data
xdata = get(dataObjs, 'XData');  
ydata = get(dataObjs, 'YData');
zdata = get(dataObjs, 'ZData');

% load data into matrix
clean_centreline = zeros(size(xdata,2),3);
clean_centreline(:,1) = xdata;
clean_centreline(:,2) = ydata;
clean_centreline(:,3) = zdata;

% Remove NaNs (created by Brushing)
clean_centreline=clean_centreline(~any(isnan(clean_centreline),2),:); 

assignin('base','cean',clean_centreline)

% Sort centreline points into order
originalx = clean_centreline(:, 1);
originaly = clean_centreline(:, 2);
originalz = clean_centreline(:, 3);
numPoints = length(originalx);
% Start with the first element (assumed start)
output = zeros(1, numPoints);
xSoFar = originalx(1);
ySoFar = originaly(1);
zSoFar = originalz(1);
remainingx = originalx(2:end); % Initialize remaining points.
remainingy = originaly(2:end); % Initialize remaining points.
remainingz = originalz(2:end); % Initialize remaining points.
for k = 1 : numPoints-1
    index = output(end);
	% Get this coordinate.
	thisX = xSoFar(k);
	thisY = ySoFar(k);
	thisZ = zSoFar(k);
	% Find distances in 3D space to all of the remaining points.
	distances = sqrt(...
		(thisX - remainingx).^2 + ...
		(thisY - remainingy).^2 + ...
		(thisZ - remainingz).^2);
	% Find the index in the remaining points of the
	% point that is closest to this one.
	[~, indexOfMin] = min(distances);
	% Assign the closest coordinates to the next point
	xSoFar(k+1) = remainingx(indexOfMin);
	ySoFar(k+1) = remainingy(indexOfMin);
	zSoFar(k+1) = remainingz(indexOfMin);
	% Remove that point from the list of remaining points.
	remainingx(indexOfMin) = [];
	remainingy(indexOfMin) = [];
	remainingz(indexOfMin) = [];
end
% Now they're arranged properly and we can proceed.

% Rename the arrays for convenience
x = xSoFar;
y = ySoFar;
z = zSoFar;

% Savitzky–Golay filter (sgolayfilt) - smoothing individual axes
windowWidth = 9;    %Standard example values
polynomialOrder = 3;
xsg=sgolayfilt(x,polynomialOrder, windowWidth);
ysg=sgolayfilt(y,polynomialOrder, windowWidth);
zsg=sgolayfilt(z,polynomialOrder, windowWidth);

% Interpolate (interparc) - standardise distance between points
t = size(xsg,2);
xyzsgpt = interparc(t, xsg, ysg, zsg);

% Plot cleaned centreline with original
hold on
plot3(x, y, z,'bx','MarkerSize', 3,'Parent',handles.clean_midline_plot);
plot3(xsg, ysg, zsg,'rx', 'MarkerSize', 5,'Parent',handles.clean_midline_plot);
plot3(xyzsgpt(:,1),xyzsgpt(:,2),xyzsgpt(:,3),'gx','Parent',handles.clean_midline_plot);
hold off
grid on
view([44,6])

sorted_centreline=zeros(size(xsg,2),3);
sorted_centreline(:,1)=xsg';
sorted_centreline(:,2)=ysg';
sorted_centreline(:,3)=zsg';

% Mean error for coordinate deviation
error_midline = (clean_centreline-sorted_centreline);%ERROR PER POINT
e_m=mean(error_midline);%ERROR PER AXIS
e_m=mean(e_m);%ERROR AVERAGE
handles.data_table(6,3)=e_m;

% Calculate parameters of centreline
handles.data_table(:,2)=centreline_tort(xyzsgpt,handles.voxel_dim);

% Error in DM
 handles.data_table(1,3) = handles.data_table(1,2)*e_m;
% Same error as distance metric for ICM
handles.data_table(3,3)= handles.data_table(1,2)*e_m;
handles.data_table(4,3)= handles.data_table(1,2)*e_m;

set(handles.tort_table,'Data',handles.data_table)

guidata(hObject, handles);


% --- Executes on slider movement.
function branch_slider_Callback(hObject, eventdata, handles)
% hObject    handle to branch_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% get branch slider values
THR = get(hObject,'Value');

% Find branch points
[~,node,link] = Skel2Graph3D(handles.skel,THR);
a=size(handles.skel,1);
b=size(handles.skel,2);
c=size(handles.skel,3);
skel_cleaned = Graph2Skel3D(node,link,a,b,c);

% Extract coordinate locations
w=size(skel_cleaned,1);
l=size(skel_cleaned,2);
h=size(skel_cleaned,3);
[x,y,z]=ind2sub([w,l,h],find(skel_cleaned(:)));
handles.x_clean = x;
handles.y_clean = y;
handles.z_clean = z;

% Plot data with removed branches
plot3(x,y,z,'square','Markersize',2,'MarkerFaceColor','b',...
    'Color','b','Parent',handles.clean_midline_plot);
view([-39,6])
fixed_xyz = [x,y,z];

% Calculate centreline parameters
handles.data_table(:,1)=centreline_tort(fixed_xyz,handles.voxel_dim);

set(handles.tort_table,'Data',handles.data_table)

guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function branch_slider_CreateFcn(hObject, eventdata, handles) %#ok<INUSD>
% hObject    handle to branch_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function centreline_button_CreateFcn(hObject, eventdata, handles)
% hObject    handle to centreline_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
