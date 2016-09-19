function [ data_array ] = centreline_tort_beta( centreline_data, voxel_dim )
% Description
% -----------
% This file contains all of the quantification methods used in the
% experiment, including the simulations
%
% Variables
% ---------
% centreline_data - centreline data [n,3] 3->[x,y,z]
% voxel_dim       - dimensions of voxel [cm]
% 
% Returns
% ---------
% data_array      - array containing all calculated variables
% data_array(1)   - distance metric
% data_array(2)   - sum of angles metric
% data_array(3)   - inflection count metric (normal)
% data_array(4)   - inflection count metric (binomial)
% data_array(5)   - length of centreline
% data_array(6)   - TBA
% data_array(7)   - TBA

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
data_xyz(:,3)=data_xyz(:,3)*voxel_dim(3);

%---------------------------------------------------------
% Length of Centreline
%---------------------------------------------------------
d = 0;
for v = 1:size(data_xyz,1)-1
    d = d + ((data_xyz(v+1,1)-data_xyz(v,1))^2 + ...
             (data_xyz(v+1,2)-data_xyz(v,2))^2 + ...
             (data_xyz(v+1,3)-data_xyz(v,3))^2 )^0.5;
end
data_table(5) = d;

%---------------------------------------------------------
%Distance Metric
%---------------------------------------------------------
start_end = [data_xyz(end-1,1)-data_xyz(1,1),...
             data_xyz(end-1,2)-data_xyz(1,2),...
             data_xyz(end-1,3)-data_xyz(1,3)];
length_start_end = (start_end(1)^2 + start_end(2)^2 + start_end(3)^2 )^0.5;
data_table(1)= d/length_start_end;

%----------------------------------------------------------
%Inflection Count Metric [t-tangent, n-normal, b-binomial]
%----------------------------------------------------------
[~,n,b]=frenet(data_xyz(:,1),data_xyz(:,2),data_xyz(:,3));
IC_n = 0;
IC_b = 0;
for i = 1:size(n,1)-2
    %Check if accellerator vector is too small [Bullitt]
    T1 = data_xyz(i+1,:) - data_xyz(i,:);
    T2 = data_xyz(i+2,:) - data_xyz(i+1,:);
    A=T2-T1;
    if norm(A)<10^-5 % in [mm]
        continue
    end
    del_n = n(i+1,:)-n(i,:);
    del_b = b(i+1,:)-b(i,:);
    local_max_n = dot(del_n,del_n);
    local_max_b = dot(del_b,del_b);
    if local_max_n > 1
        IC_n = IC_n + 1;
    end
    if local_max_b > 1
        IC_b = IC_b + 1;
    end
end
data_table(3)= (IC_n+1)*(d/length_start_end)/10;
data_table(4)= (IC_b+1)*(d/length_start_end)/10;

%---------------------------------------------------------------
% Sum Of Angles Metric [SOAM]
%---------------------------------------------------------------
n = size(data_xyz,1);
CP = zeros(1,n);
% Consecutive 3 points on centreline
for k = 2:n-3
    % Calculate total angle from in-plane angle and torsion angle
    T1 = data_xyz(k,:) - data_xyz(k-1,:);
    T2 = data_xyz(k+1,:) - data_xyz(k,:);
    T3 = data_xyz(k+2,:) - data_xyz(k+1,:);
    
    % In-plane angle
    IPA = acos(dot(T1/norm(T1),T2/norm(T2)));
    
    % Torsion angle
    cross_mag_T12 = cross(T1,T2);
    cross_mag_T12 = sqrt(cross_mag_T12(1)^2 + cross_mag_T12(2)^2 + cross_mag_T12(3)^2);
    cross_mag_T23 = cross(T2,T3);
    cross_mag_T23 = sqrt(cross_mag_T23(1)^2 + cross_mag_T23(2)^2 + cross_mag_T23(3)^2);
    a=abs(cross(T1,T2))/abs(cross_mag_T12);
    b=abs(cross(T2,T3))/abs(cross_mag_T23);
    TPA = acos(dot(a,b));  %#ok<NASGU>
    TPA = 0; % Ignore torsion angle
    
    % Total angle
    CP(k) = sqrt((IPA^2)+(TPA^2));
end

% Remove NaNs
CP_NaN = CP;
CP_NaN(isnan(CP_NaN)) = 0;
fixed_NaN = data_xyz;

% Sum of Angles Metric Equation
n = size(fixed_NaN,1);
a = CP_NaN(1:(n-3));
b = norm(fixed_NaN(2:n,:)-fixed_NaN(1:n-1,:));
SOAM = sum(a)/(sum(b)); 
data_table(2) = real(SOAM)/1000;

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

data_table(7) = norm(k)/1000;

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

