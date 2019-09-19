%% Assignment 2 Tom Manuel
close all
clear all
load('pointTargetData.mat')
data = veraStrct.data;
%imagesc(data(:,:,64),[-100 100])

%% Part 4: single fixed focus beamforming
%pixel in physical space
close all
Nx = veraStrct.numElementsPerXmt;
dx = 1E-3 * veraStrct.XMTspacingMM;

%get Nx,dx,Nz, and dz
foo = size(data);
t0 = veraStrct.timeZero -1; % nPts to throw away 
Nz = foo(1)-t0;
data = data(t0+1:end,:,:);
% woah do you have to guess c to get dz??
% dz = [m/s]*[s/sample] = c * fs^-1 
c = 1500; %m/s
dz = .5 * c / 20E6;

% build coordinate matrices
Xe = dx.*linspace(-Nx/2,Nx/2,Nx);
Ze = 0;
% set focus (in meters)
xf = 0;
zf = .04; %4 cm focal depth

% calculate delay
Trx = (1/c)*sqrt((Xe-xf).^2 + (Ze-zf).^2) + zf/c; %[s]
% subtract max
Trx = Trx -max(Trx);

% apply the delay to a single image
% make a time vector to use with interp1
fs = veraStrct.samplingRateMHz;
t = 1E-6 .*linspace(0,Nz/fs,Nz);
im0 = data(:,:,64); % original data for one frame
imd = zeros(size(im0)); % delayed frame
for i=1:128
    imd(:,i) = interp1(t,im0(:,i),t+Trx(i),'linear');
end

% apply the delay to all frames
%exploit the fact that for all beams, the distribution of delay is the same
datad = zeros(size(data));
for i=1:128
    datad(:,i,:) = interp1(t,data(:,i,:),t+Trx(i),'linear');
end

Zvec = dz.*linspace(0,Nz-1,Nz);
% sum
imf = squeeze(sum(datad,2));
imf(isnan(imf))=-100;
% compress
imf = 20.*log10(abs(hilbert(imf)));

%% Part 4 deliverables Point target
% image of original channel data
% make vector representing depth

figure
subplot(1,3,1)
imagesc(Xe.*1000,Zvec.*1000,im0(500:1200,:),[10 100])
colormap('gray')
title('Original channel data (zoomed and saturated)')
xlabel('Xe (mm)')
ylabel('Depth (mm)')

% image of delayed channel data
subplot(1,3,2)
imagesc(Xe.*1000,Zvec.*1000,imd(500:1200,:),[10 100])
colormap('gray')
title('Delayed channel data (zoomed and saturated)')
xlabel('Xe (mm)')
ylabel('Depth (mm)')

% image of delayed and summed data
subplot(1,3,3)
imagesc(Xe.*1000,Zvec.*1000,imf,[10 100])
colormap('gray')
title('Delayed and Summed')
xlabel('Xe (mm)')
ylabel('Depth (mm)')


%% Repeat all steps but for cyst
close all
clear all
load('anecoicCystData.mat')
data = veraStrct.data;

close all
Nx = veraStrct.numElementsPerXmt;
dx = 1E-3 * veraStrct.XMTspacingMM;

%get Nx,dx,Nz, and dz
foo = size(data);
t0 = veraStrct.timeZero -1; % nPts to throw away 
Nz = foo(1)-t0;
data = data(t0+1:end,:,:);
% woah do you have to guess c to get dz??
% dz = [m/s]*[s/sample] = c * fs^-1 
c = 1500; %m/s
dz = .5 * c / 20E6;

% build coordinate matrices
Xe = dx.*linspace(-Nx/2,Nx/2,Nx);
Ze = 0;
% set focus (in meters)
xf = 0;
zf = .04; %4 cm focal depth

% calculate delay
Trx = (1/c)*sqrt((Xe-xf).^2 + (Ze-zf).^2) + zf/c; %[s]
% subtract max
Trx = Trx -max(Trx);

% apply the delay to a single image
% make a time vector to use with interp1
fs = veraStrct.samplingRateMHz;
t = 1E-6 .*linspace(0,Nz/fs,Nz);
im0 = data(:,:,64); % original data for one frame
imd = zeros(size(im0)); % delayed frame

%%
for i=1:128
    %shift with interp 1
    %imd(:,i) = interp1(t,im0(:,i),t+Trx(i),'linear');
    
    %shift discretely
    
end

% apply the delay to all frames
%exploit the fact that for all beams, the distribution of delay is the same
datad = zeros(size(data));
for i=1:128
    datad(:,i,:) = interp1(t,data(:,i,:),t+Trx(i),'linear');
end

Zvec = dz.*linspace(0,Nz-1,Nz);
% sum
imf = squeeze(sum(datad,2));
imf(isnan(imf))=-100;
% compress
imf = 20.*log10(abs(hilbert(imf)));



%% Part 4 deliverables Cyst
% image of original channel data
% make vector representing depth

figure
% image of delayed and summed data
imagesc(Xe.*1000,Zvec.*1000,imf,[10 100])
colormap('gray')
title('Delayed and summed')
xlabel('Xe (mm)')
ylabel('Depth (mm)')
axis image


%%
% may be useful for adaptive focusing

% %% Apply delay to a single image
% 
% %First I will create a grid that represents the x and z coordinates of each
% %pixel in physical space
% Nx = veraStrct.numElementsPerXmt;
% dx = veraStrct.XMTspacingMM;
% 
% %get Nx,dx,Nz, and dz
% foo = size(data);
% Nz = foo(1);
% % woah do you have to guess c to get dz??
% % dz = [m/s]*[s/sample] = c * fs^-1 
% c = 1500; %m/s
% dz = c / 20E6;
% 
% % build coordinate matrices
% Xe = repmat(dx.*linspace(-Nx/2,Nx/2,Nx),[Nz 1]);
% Ze = repmat(dz.*linspace(0,Nz-1,Nz)',[1,Nx]);
% 
% % set focus (in meters)
% xf = 0;
% zf = .04;
% 
% % calculate delay
% Trx = (1/c)*sqrt((Xe-xf).^2 + (Ze-zf).^2) + zf/c;





