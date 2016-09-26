function [ data_array ] = centreline_tort_beta( centreline_data, voxel_dim, spacing )
% Description
% -----------
% This file contains all of the quantification methods used in the
% experiment, including the simulations
%
% Variables
% ---------
% centreline_data - centreline data [n,3] 3->[x,y,z]
% voxel_dim       - [x,y,z] dimensions of voxel [metres]
% spacing         - spacing between slices [metres]
% Returns
% ---------
% data_array      - array containing all calculated variables
% data_array(1)   - distance metric
% data_array(2)   - sum of angles metric
% data_array(3)   - inflection count metric (normal)
% data_array(4)   - inflection count metric (binomial)
% data_array(5)   - length of centreline
% data_array(6)   - standard deviation of average curvature
% data_array(7)   - norm of curvature

%---------------------------------------------------------
% Initial Processing
%---------------------------------------------------------
% initialize input data array
data_xyz=centreline_data;

% initialize output array
data_table = zeros(1,7);

% multiply centreline dimensions by voxel dimensions
data_xyz(:,1)=data_xyz(:,1)*voxel_dim(1);
data_xyz(:,2)=data_xyz(:,2)*voxel_dim(2);
data_xyz(:,3)=data_xyz(:,3)*(voxel_dim(3)+spacing);

%---------------------------------------------------------
% Length of Centreline [mm]
%---------------------------------------------------------
d = 0;
for v = 1:size(data_xyz,1)-1
    d = d + ((data_xyz(v+1,1)-data_xyz(v,1))^2 + ...
             (data_xyz(v+1,2)-data_xyz(v,2))^2 + ...
             (data_xyz(v+1,3)-data_xyz(v,3))^2 )^0.5;
end
data_table(5) = d*10^3;

%---------------------------------------------------------
% Distance Metric
%---------------------------------------------------------
start_end = [data_xyz(end-1,1)-data_xyz(1,1),...
             data_xyz(end-1,2)-data_xyz(1,2),...
             data_xyz(end-1,3)-data_xyz(1,3)];
length_start_end = (start_end(1)^2 + start_end(2)^2 + start_end(3)^2 )^0.5;
data_table(1)= d/length_start_end;

%----------------------------------------------------------
% Inflection Count Metric [t-tangent, n-normal, b-binomial]
%----------------------------------------------------------
P=data_xyz;
% Calculate frenet frame of centreline
[~,n,b]=frenet(data_xyz(:,1),data_xyz(:,2),data_xyz(:,3));
IC_n = 0;
IC_b = 0;
for i = 2:size(n,1)-1
    
    % First check if accellerator vector is too small [Bullitt]
    T1 = P(i,:) - P(i-1,:);
    T2 = P(i+1,:) - P(i,:);
    A=T2-T1; % Acceleration vector
    if norm(A)<10^-5 % in [mm]???
        continue
    end
    
    % Calculate inflection count (del_ : change in axis)
    del_n = n(i,:)-n(i-1,:);
    del_b = b(i,:)-b(i-1,:);
    local_max_n = dot(del_n,del_n);
    local_max_b = dot(del_b,del_b);
    
    if local_max_n > 1 % normal axis maxima
        IC_n = IC_n + 1;
    end
    if local_max_b > 1 %binomial axis maxima
        IC_b = IC_b + 1;
    end
end
data_table(3)= (IC_n+1)*(d/length_start_end);
data_table(4)= (IC_b+1)*(d/length_start_end);

%---------------------------------------------------------------
% Sum Of Angles Metric [SOAM] UNITS: radians/cm
%---------------------------------------------------------------
n = size(data_xyz,1);
CP = zeros(1,n);
P = data_xyz; % Convert into cm
% Consecutive 3 points on centreline
for k = 2:n-3
    % Calculate total angle from in-plane angle and torsion angle
    
    % Initial difference vectors
    T1 = P(k,:) - P(k-1,:);
    T2 = P(k+1,:) - P(k,:);
    
    % In-plane angle
    IP = acos(dot(T1/norm(T1),T2/norm(T2)));
    
    % Torsion angle ignored as recommended by Bullitt
    % Only introduced noise in previous studies
    
    % Total angle
    CP(k) = sqrt((IP^2));
end

% Remove NaNs
CP(isnan(CP)) = 0;

% Sum of Angles Metric Equation
n = size(P,1);
SOAM = sum(CP(1:(n-3)))/(sum(norm(P(2:n,:)-P(1:n-1,:)))); 
data_table(2) = real(SOAM)*(pi/180); %convert into radians / mm

%---------------------------------------------------------------
% Standard deviation of curvature [SDOC]
%---------------------------------------------------------------

% smoothing parameter
dS = 0.5*10^-6;

% extract frenet frame - tangent
[~,~,~,k,~]=frenet(data_xyz(:,1),data_xyz(:,2),data_xyz(:,3));

% curvature average length
k_m = (dS/d)*sum(k(3:size(data_xyz,1)-2));

% curvature standard deviation
k_sd = ((dS/d)^0.5) * sqrt(sum((k(3:size(data_xyz,1)-2)-k_m).^2));

data_table(6)= k_sd;

%------------------------------------------------------------
% Norm of Curvature
%------------------------------------------------------------
dS = 0.262*10^-3; % played with parameter
data_table(7) = norm(k,1)*dS;

% Load data into array
data_array=data_table;


% Frenet frame as the basis of curvature information
% How to use this to get reliables frequency info...
% Try t b n each and see across different parameters
%


% Volume of Tube Segment T.B.C.
% faces = get(handles.vesselpatch, 'Faces'); 
% vertices = get(handles.vesselpatch, 'Vertices'); 
% volume = meshVolume(vertices,faces);
% volume = abs(volume*handles.voxel_dim(1)*handles.voxel_dim(2)*(handles.voxel_dim(3)+handles.spacing));
% handles.data_table(7,1) = volume;
end

