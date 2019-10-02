%% Tom
% tired of not understanding these examples
% im starting from scratch here

% make a single point target
amp = 1;
p = [0 0 50]./1000; % 70mm

% build a transducer to transmit and receive

% Set initial parameters
f0=3e6; % Transducer center frequedncy [Hz]
fs=100e6; % Sampling frequency [Hz]
c=1540; % Speed of sound [m/s]
lambda=c/f0; % Wavelength [m]
height=5/1000; % Height of element [m]
width=lambda; % Width of element [m]
kerf=width/100; % Distance between transducer elements [m]
N_elements=128; % Number of elements

focus=[0 0 70]/1000; % Initial electronic focus

set_field('c',c);
set_field('fs',fs);

% Define the transducers
TE = xdc_linear_array (N_elements, width, height, kerf, 1, 1, focus); %emit transducer
TR = xdc_linear_array (N_elements, width, height, kerf, 1, 1, focus); %receive transducer
% Set the impulse response and excitation of the emit aperture

% set minimimum quantization
min_delay = (1/f0)/20;
xdc_quantization (TE, min_delay);
xdc_quantization (TR, min_delay);

impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse (TE, impulse_response); %emit
xdc_impulse (TR, impulse_response); %receive
excitation=sin(2*pi*f0*(0:1/fs:2/f0));
excitation = [excitation,0];
xdc_excitation (TE, excitation);

xdc_focus (TR, 0, p);

% scan a the sample
no_lines = 20;
image_width=no_lines/1000;  
d_x = image_width/no_lines;
image_data=0;
x = -image_width/2;
p(:,1) = x;

for i=1:no_lines
    [v,t]=calc_scat(TE, TR, p, amp);
    image_data(1:max(size(v)),i)=v;
    times(i) = t;
    
    % slide p to the left
    p(:,1) = p(:,1) + d_x;
end

mk_img_tom

%% plot lateral psf
% image is in new_env

% get row with max
[~,I] = max(new_env(:,50));

hold off
lpsf = new_env(I,:);
plot(xx,lpsf)


%% we're good now actually do the project
foo = 511; %helps make figure
close all
f0=3e6;
Quants = (1/f0)./[10 5 1.5 1.2 1];
f1 = figure;
f2 = figure;
for q = 1:length(Quants)
    
    % make a single point target
    amp = 1;
    p = [0 0 50]./1000; % 70mm

    % build a transducer to transmit and receive

    % Set initial parameters
    f0=3e6; % Transducer center frequedncy [Hz]
    fs=100e6; % Sampling frequency [Hz]
    c=1540; % Speed of sound [m/s]
    lambda=c/f0; % Wavelength [m]
    height=5/1000; % Height of element [m]
    width=lambda; % Width of element [m]
    kerf=width/100; % Distance between transducer elements [m]
    N_elements=128; % Number of elements

    focus=[0 0 70]/1000; % Initial electronic focus

    set_field('c',c);
    set_field('fs',fs);

    % Define the transducers
    TE = xdc_linear_array (N_elements, width, height, kerf, 1, 1, focus); %emit transducer
    TR = xdc_linear_array (N_elements, width, height, kerf, 1, 1, focus); %receive transducer
    % Set the impulse response and excitation of the emit aperture

    % set minimimum quantization
    min_delay = Quants(q);
    xdc_quantization (TE, min_delay);
    xdc_quantization (TR, min_delay);

    impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
    impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
    xdc_impulse (TE, impulse_response); %emit
    xdc_impulse (TR, impulse_response); %receive
    excitation=sin(2*pi*f0*(0:1/fs:2/f0));
    excitation = [excitation,0];
    xdc_excitation (TE, excitation);

    xdc_focus (TR, 0, p);

    % scan a the sample
    no_lines = 20;
    image_width=no_lines/1000;  
    d_x = image_width/no_lines;
    image_data=0;
    x = -image_width/2;
    p(:,1) = x;

    for i=1:no_lines
        [v,t]=calc_scat(TE, TR, p, amp);
        image_data(1:max(size(v)),i)=v;
        times(i) = t;

        % slide p to the left
        p(:,1) = p(:,1) + d_x;
    end

    mk_img_tom
    
    % show psfs
    figure(f1)
    subplot(foo)
    imagesc(xx,zz(150:180),new_env(150:180,:),[-50 0])
    axis image
    colormap(gray)
    foo = foo +1;
    axis off
    set(gcf,'color','black');
    
    %plot psfs
    %get row with max
    figure(f2)
    [~,I] = max(new_env(:,50));
    lpsf = new_env(I,:);
    hold on
    plot(xx,lpsf,'linewidth',1.5)
    
end

%%
figure(f2)
legend('\lambda/10','\lambda/5', '\lambda/1.5', '\lambda/1.2','\lambda')
title('Lateral PSFs with increasing quantization')
xlabel('X (mm)')
ylabel('dB')





