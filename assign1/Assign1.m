% Assignment 1 Problem 4 & 5
% Tom Manuel

%% Part A
f0 = 1E6;
fs = 40E6;
c=1540;

%create time
%first figure out how many points to use in time as to cover
%exactly 3 cycles
ncyc = 3;
T = 1/f0;
npts = ncyc*T*fs; %number of points
t = linspace(0,T*3,npts);

%point along x axis (m)
x = 2E-6;

%get wavenumber (k)
%k = 2pi /lamda
k = 2 * pi / (c/f0 );

%create sinusoid
% A*sin(wt+kx)
X = 1*sin(2*pi*f0*t+ k*x);
close all
plot(t,X)
xlabel('time (s)')
ylabel('amplitude')

%% Part B

%sig from Part A
X = 1*sin(2*pi*f0*t+ k*x);

% create same signal shifted Pi/8
phase = pi/8;
X2 = 1*sin(2*pi*f0*t+ k*x + phase);
close all
plot(t,X,t,X2)
legend('X','X shifted')
xlabel('time (s)')
ylabel('amplitude')

%% Part C

% get the time value corresponded to a phase shift of pi/8
% 2pi = 1/f0, 2pi/pi/8= 16, shift = (1/f0)/16
shift = -(1/f0)/16;

Xnear = interp1(t,X2,t+shift,'nearest');
Xlin = interp1(t,X2,t+shift,'linear');
Xspline = interp1(t,X2,t+shift,'spline');
close all
plot(t,X,t,X2,t,Xnear,t,Xlin,t,Xspline)
legend('X','X2','Nearest','Linear','Spline')
xlim([2E-6 2.5E-6])
xlabel('time (s)')
ylabel('amplitude')

%% Part D

% Compute mean square error 
errN = nanmean((X-Xnear).^2) %nearest
errL = nanmean((X-Xlin).^2) %linear
errS = nanmean((X-Xspline).^2) %spline

%calculate power in signal
%amplitude was 1
p = 1/sqrt(2);
dBN = 10*log10(errN/p)
dBL = 10*log10(errL/p)
dBS = 10*log10(errS/p)

%% 5a Ex:6.1(page 62 USERGUIDE)

% Example of use of the new Field II program running under Matlab
% This example shows how a phased array B-mode system scans an image
% This script assumes that the field_init procedure has been called
% Example by Joergen Arendt Jensen, Nov. 28, 1995.

% Generate the transducer apertures for send and receive
f0=3e6; % Transducer center frequency [Hz]
fs=100e6; % Sampling frequency [Hz]
c=1540; % Speed of sound [m/s]
lambda=c/f0; % Wavelength
element_height=5/1000; % Height of element [m]
kerf=0.1/1000; % Kerf [m]
focus=[0 0 70]/1000; % Fixed focal point [m]

% Generate aperture for emission
emit_aperture = xdc_linear_array (128, lambda/2, element_height, kerf, 1, 1,focus);

% Set the impulse response and excitation of the emit aperture
impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse (emit_aperture, impulse_response);
excitation=sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation (emit_aperture, excitation);

% Generate aperture for reception
receive_aperture = xdc_linear_array (128, lambda/2, element_height, kerf, 1, 1,focus);

% Set the impulse response for the receive aperture
xdc_impulse (receive_aperture, impulse_response);

% Do phased array imaging
point_position=[0 0 70]/1000; % Position of the point to be imaged
no_lines=50; % Number of A-lines in image
sector=20 * pi/180; % Size of image sector
d_theta=sector/no_lines; % Increment in angle for 90 deg. image

% Pre-allocate some storage
image_data=zeros(800,no_lines);
theta= -sector/2;

for i=1:no_lines
    % Set the focus for this direction
    xdc_focus (emit_aperture, 0, [70*sin(theta) 0 70*cos(theta)]/1000);
    xdc_focus (receive_aperture, 0, [70*sin(theta) 0 70*cos(theta)]/1000);
    % Calculate the received response
    [v, t1]=calc_scat(emit_aperture, receive_aperture, point_position, 1);
    % Store the result
    image_data(1:max(size(v)),i)=v';
    times(i) = t1;
    % Steer in another angle
    theta = theta + d_theta;
end
% Here the display of the data is inserted
plot(image_data)

%% 5b Ex:6.2

% Example of use of the new Field II program running under Matlab
% This example shows how a linear array B-mode system scans an image
% This script assumes that the field_init procedure has been called
% Example by Joergen Arendt Jensen, Version 2.0, March 22, 2011.
% Generate the transducer apertures for send and receive

f0=3e6; % Transducer center frequency [Hz]
fs=100e6; % Sampling frequency [Hz]
c=1540; % Speed of sound [m/s]
lambda=c/f0; % Wave length [m]
width=lambda; % Width of element
element_height=5/1000; % Height of element [m]
kerf=width/20; % Kerf [m]
focus=[0 0 50]/1000; % Fixed focal point [m]
N_elements=192; % Number of elements in the transducer
N_active=64; % Active elements in the transducer
% Set the sampling frequency
set_sampling(fs);
% Generate aperture for emission
emit_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 5, focus);
% Set the impulse response and excitation of the emit aperture
impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse (emit_aperture, impulse_response);
excitation=sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation (emit_aperture, excitation);
% Generate aperture for reception
receive_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 5, focus);
% Set the impulse response for the receive aperture
xdc_impulse (receive_aperture, impulse_response);
% Load the computer phantom
[phantom_positions, phantom_amplitudes] = cyst_phantom(10000);
% Do linear array imaging
no_lines=N_elements-N_active+1; % Number of A-lines in image
dx=width; % Increment for image
z_focus=50/1000;
% Pre-allocate some storage
image_data=zeros(1,no_lines);
for i=1:no_lines
    % Find position for imaging
    x=(i-1-no_lines/2)*dx;
    % Set the focus for this direction
    xdc_center_focus (emit_aperture, [x 0 0]);
    xdc_focus (emit_aperture, 0, [x 0 z_focus]);
    xdc_center_focus (receive_aperture, [x 0 0]);
    xdc_focus (receive_aperture, 0, [x 0 z_focus]);
    % Set the active elements using the apodization
    apo=[zeros(1, i-1) hamming(N_active)' zeros(1, N_elements-N_active-i+1)];
    xdc_apodization (emit_aperture, 0, apo);
    xdc_apodization (receive_aperture, 0, apo);
    % Calculate the received response
    [v, t1]=calc_scat(emit_aperture, receive_aperture, phantom_positions, phantom_amplitudes);
    % Store the result
    image_data(1:max(size(v)),i)=v;
    times(i) = t1;
end
% Free space for apertures
xdc_free (emit_aperture)
xdc_free (receive_aperture)
% Adjust the data in time and display it as
% a gray scale image
min_sample=min(times)*fs;
for i=1:no_lines
    rf_env=abs(hilbert([zeros(round(times(i)*fs-min_sample),1); image_data(:,i)]));
    env(1:size(rf_env,1),i)=rf_env;
end

% make logarithmic compression to a 60 dB dynamic range
% with proper units on the axis
env_dB=20*log10(env);
env_dB=env_dB-max(max(env_dB));
env_gray=127*(env_dB+60)/60;
depth=((0:size(env,1)-1)+min_sample)/fs*c/2;
x=((1:no_lines)-no_lines/2)*dx;
image(x*1000, depth*1000, env_gray)
xlabel('Lateral distance [mm]')
ylabel('Depth [mm]')
axis('image')
colormap(gray(128))
title('Image of cyst phantom (60 dB dynamic range)')
