%% Assignment 2 Tom Manuel
% Part 4
close all
clear all
load('pointTargetData.mat')
data = veraStrct.data;
%imagesc(data(:,:,64),[-100 100])

% Part 4: single fixed focus beamforming
%pixel in physical space
close all
Nx = veraStrct.numElementsPerXmt;
dx = 1E-3 * veraStrct.XMTspacingMM;

%get Nx,dx,Nz, and dz
foo = size(data);
t0 = veraStrct.timeZero -1; % nPts to throw away 
Nz = foo(1)-t0;
data = data(t0+1:end,:,:);
foo = size(data);
% woah do you have to guess c to get dz in [m]??
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
Trx = Trx - min(Trx); %[s] add 1 to take away that pesky zero
% since shifting discretely, upsample data
usf = 8;
z1 = linspace(1,foo(1),foo(1));
z2 = linspace(1,foo(1),foo(1)*usf);
datau = interp1(z1,data,z2,'linear'); %upsampled data

%convert Trx from [s] to samples
%trx = trx [s] * fs [sam/s] *usf [sam/sam]
%round to descritize
fs = veraStrct.samplingRateMHz * 1E6;
Trx = round(Trx.*fs.*usf) +1; %add 1 to get rid of pesky zero
%
% apply the delay to a single image
im0 = datau(:,:,64); % original data for one frame
imd = zeros(size(im0)); % delayed frame
for i=1:128
    imd(1:end-Trx(i)+1,i) = im0(Trx(i):end,i);
end

% apply the delay to all frames
%exploit the fact that for all beams, the distribution of delay is the same
datad = zeros(size(datau)); %data delayed
for i=1:128
    %datad(:,i,:) = interp1(t,data(:,i,:),t+Trx(i),'linear');
    datad(1:end-Trx(i)+1,i,:) = datau(Trx(i):end,i,:);
end


Zvec = dz.*linspace(0,Nz-1,Nz.*usf); % z axis for plotting time
% sum
imf = squeeze(sum(datad,2));
imf(isnan(imf))=-100;
% compress
imf = 20.*log10(abs(hilbert(imf)));

% Part 4 deliverables Point target
% image of original channel data
% make vector representing depth

figure
subplot(1,3,1)
imagesc(im0(900*usf:1100*usf,:),[10 100])
colormap('gray')
title('Original channel data (zoomed and saturated)')
xlabel('Xe (px)')
ylabel('Depth (px)')

% image of delayed channel data
subplot(1,3,2)
imagesc(imd(900*usf:1100*usf,:),[10 100])
colormap('gray')
title('Delayed channel data (zoomed and saturated)')
xlabel('Xe (px)')
ylabel('Depth (px)')

% image of delayed and summed data
subplot(1,3,3)
imagesc(Xe.*1000,Zvec.*1000,imf,[10 100])
colormap('gray')
title('Delayed and Summed')
xlabel('Xe (mm)')
ylabel('Depth (mm)')
axis image


%% Still Part 4
% Repeat all steps but for cyst
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
foo = size(data);
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
Trx = Trx - min(Trx); %[s] add 1 to take away that pesky zero
% since shifting discretely, upsample data
usf = 10;
z1 = linspace(1,foo(1),foo(1));
z2 = linspace(1,foo(1),foo(1)*usf);
datau = interp1(z1,data,z2,'linear'); %upsampled data

%convert Trx from [s] to samples
%trx = trx [s] * fs [sam/s] *usf [sam/sam]
%round to descritize
fs = veraStrct.samplingRateMHz * 1E6;
Trx = round(Trx.*fs.*usf) +1; %add 1 to get rid of pesky zero

% apply the delay to all frames
%exploit the fact that for all beams, the distribution of delay is the same
datad = zeros(size(datau)); %data delayed
for i=1:128
    %datad(:,i,:) = interp1(t,data(:,i,:),t+Trx(i),'linear');
    datad(1:end-Trx(i)+1,i,:) = datau(Trx(i):end,i,:);
end

Zvec = dz.*linspace(0,Nz-1,Nz.*usf); % z axis for plotting time
% sum
imf = squeeze(sum(datad,2));
imf(isnan(imf))=-100;
% compress
imf = 20.*log10(abs(hilbert(imf)));

% Part 4 deliverables Cyst
% image of original channel data
% make vector representing depth

figure
% image of delayed and summed data
Zvec = dz.*linspace(0,Nz-1,Nz); % z axis for plotting time
imagesc(Xe.*1000,Zvec.*1000,imf,[10 100])
colormap('gray')
title('Delayed and summed')
xlabel('Xe (mm)')
ylabel('Depth (mm)')
axis image


%% Part 5, Continuous delay and sum
close all
clear all
load('pointTargetData.mat')
data = veraStrct.data;
%imagesc(data(:,:,64),[-100 100])

%First I will create a grid that represents the x and z coordinates of each
%pixel in physical space
Nx = veraStrct.numElementsPerXmt;
dx = 1E-3*veraStrct.XMTspacingMM;

%get Nx,dx,Nz, and dz
foo = size(data);
Nz = foo(1);
t0 = veraStrct.timeZero -1; % nPts to throw away 
Nz = foo(1)-t0;
data = data(t0+1:end,:,:);
foo = size(data);
% woah do you have to guess c to get dz??
% dz = [m/s]*[s/sample] = c * fs^-1 
c = 1500; %m/s
dz = .5 * c / 20E6;

% build coordinate matrices
Xe = repmat(dx.*linspace(-Nx/2,Nx/2,Nx),[Nz 1]);
Ze = repmat(dz.*linspace(0,Nz-1,Nz)',[1,Nx]);

% set focus (in meters)
xf = 0;

% create a time vector
fs = veraStrct.samplingRateMHz;
t = 1E-6 .*linspace(0,Nz/fs,Nz);

% calculate delay
Trx = (1/c)*sqrt((Xe-xf).^2 + (Ze).^2) + Ze/c; %[s]
%Trx = Trx - min(Trx(:));
Trx = Trx - repmat(min(Trx,[],2),[1 Nx]);
%
h=waitbar(0,'processing');
datad = zeros(size(data)); %delayed data
for i=1:128
    datad(:,i,:) = interp1(t,data(:,i,:),t+Trx(:,i)','linear');
    waitbar(i/128,h);
end 
close(h);

% sum
imf = squeeze(sum(datad,2)); %sum
imf(isnan(imf))=-1000; %interp1 drops some Nans in here
% compress
imf = 20.*log10(abs(hilbert(imf)));

% image of delayed and summed data
Zvec = dz.*linspace(0,Nz-1,Nz); % z axis for plotting time
figure
imagesc(Xe(1,:).*1000,Zvec.*1000,imf,[10 100])
colormap('gray')
title('Dynamically focused')
xlabel('Xe (mm)')
ylabel('Depth (mm)')
axis image

%% Part 5b, Continuous delay and sum (cyst data)
close all
clear all
load('anecoicCystData.mat')
data = veraStrct.data;
%imagesc(data(:,:,64),[-100 100])

%First I will create a grid that represents the x and z coordinates of each
%pixel in physical space
Nx = veraStrct.numElementsPerXmt;
dx = 1E-3*veraStrct.XMTspacingMM;

%get Nx,dx,Nz, and dz
foo = size(data);
Nz = foo(1);
t0 = veraStrct.timeZero -1; % nPts to throw away 
Nz = foo(1)-t0;
data = data(t0+1:end,:,:);
foo = size(data);

% woah do you have to guess c to get dz??
% dz = [m/s]*[s/sample] = c * fs^-1 
c = 1500; %m/s
dz = .5 * c / 20E6;

% build coordinate matrices
Xe = repmat(dx.*linspace(-Nx/2,Nx/2,Nx),[Nz 1]);
Ze = repmat(dz.*linspace(0,Nz-1,Nz)',[1,Nx]);

% set focus (in meters)
xf = 0;

% create a time vector
fs = veraStrct.samplingRateMHz;
t = 1E-6 .*linspace(0,Nz/fs,Nz);

% calculate delay
Trx = (1/c)*sqrt((Xe-xf).^2 + (Ze).^2) + Ze/c; %[s]
%Trx = Trx - min(Trx(:));
Trx = Trx - repmat(min(Trx,[],2),[1 Nx]);
%
h=waitbar(0,'processing');
datad = zeros(size(data)); %delayed data
for i=1:128
    datad(:,i,:) = interp1(t,data(:,i,:),t+Trx(:,i)','linear');
    waitbar(i/128,h);
end 
close(h);

% sum
imf = squeeze(sum(datad,2)); %sum
imf(isnan(imf))=-1000; %interp1 drops some Nans in here
% compress
imf = 20.*log10(abs(hilbert(imf)));

% image of delayed and summed data
Zvec = dz.*linspace(0,Nz-1,Nz); % z axis for plotting time
figure
imagesc(Xe(1,:).*1000,Zvec.*1000,imf,[10 100])
colormap('gray')
title('Dynamically focused')
xlabel('Xe (mm)')
ylabel('Depth (mm)')
axis image

%% Part 6 Parallel imaging
% adapted from dynamic focus code... as if that weren't complicated enough
% each foci will have a slightly different delay profile
% so calculate a Trx for each
close all
clear all
load('pointTargetData.mat')

pff = 16; % parallel focus factor (2,4,8,16)
data = veraStrct.data;
%First I will create a grid that represents the x and z coordinates of each
%pixel in physical space
Nx = veraStrct.numElementsPerXmt;
dx = 1E-3*veraStrct.XMTspacingMM;

%get Nx,dx,Nz, and dz
foo = size(data);
Nz = foo(1);
t0 = veraStrct.timeZero -1; % nPts to throw away 
Nz = foo(1)-t0;
data = data(t0+1:end,:,:);
foo = size(data);
% woah do you have to guess c to get dz??
% dz = [m/s]*[s/sample] = c * fs^-1 
c = 1500; %m/s
dz = .5 * c / 20E6;

% build coordinate matrices
Xe = repmat(dx.*linspace(-Nx/2,Nx/2,Nx),[Nz 1]);
Ze = repmat(dz.*linspace(0,Nz-1,Nz)',[1,Nx]);

% create lateral foci vector
Xf = linspace(-dx*pff/2,dx*pff/2,pff);

Trx = zeros([Nz,Nx,pff]);
% calculate delay for each focus
for i=1:pff
    ttemp = (1/c)*sqrt((Xe-Xf(i)).^2 + (Ze).^2) + Ze/c; %[s] temp delay matrix [Nz x Nx]
    %Trx = Trx - min(Trx(:));
    Trx(:,:,i) = ttemp - repmat(min(ttemp,[],2),[1 Nx]);
end

% create a time vector
fs = veraStrct.samplingRateMHz;
t = 1E-6 .*linspace(0,Nz/fs,Nz);

% takeaway beams to artificially simulate parallel foci
nBeams = foo(3);
datar = zeros([Nz,Nx,nBeams/pff]); %reduced data [Nz,Nx,nbeams/par_foc_factor]
f=1;
for i=1:pff:nBeams
    datar(:,:,f)=data(:,:,i);
    f=f+1;
end

%Now loop through datar and apply delays
bar = nBeams/pff;
datad = zeros([Nz Nx bar pff]); %[Nz Nx nBeams/pff pff]

%
bar = nBeams/pff;
for i=1:nBeams
    for j=1:pff
        datad(:,i,:,j) = interp1(t,datar(:,i,:),t+squeeze(Trx(:,i,j))','linear');
        %datad(:,i,:,j) = datar(:,i,:);
    end
end

datas = squeeze(sum(datad,2));

foo = size(data);
imf = zeros([foo(1) foo(3)]);%original data size without element dimension
for i=1:bar
    for j=1:pff
        imf(:,i*pff -pff+j) = datas(:,i,j);
    end
end

imf(isnan(imf))=-1000; %interp1 drops some Nans in here
% compress
imf = 20.*log10(abs(hilbert(imf)));

% image of delayed and summed data
Zvec = dz.*linspace(0,Nz-1,Nz); % z axis for plotting time
figure
imagesc(Xe(1,:).*1000,Zvec.*1000,imf,[10 100])
colormap('gray')
title('Parallel Imaging')
xlabel('Xe (mm)')
ylabel('Depth (mm)')
axis image

%% display 4 images with different pff
%pIms = zeros([2353,128,4]);
figure
subplot(1,4,1)
imagesc(Xe(1,:).*1000,Zvec.*1000,pIms(:,:,1),[10 100])
colormap('gray')
title('2 parallel foci')
xlabel('Xe (mm)')
ylabel('Depth (mm)')
axis image

subplot(1,4,2)
imagesc(Xe(1,:).*1000,Zvec.*1000,pIms(:,:,2),[10 100])
colormap('gray')
title('4 parallel foci')
xlabel('Xe (mm)')
ylabel('Depth (mm)')
axis image

subplot(1,4,3)
imagesc(Xe(1,:).*1000,Zvec.*1000,pIms(:,:,3),[10 100])
colormap('gray')
title('8 parallel foci')
xlabel('Xe (mm)')
ylabel('Depth (mm)')
axis image

subplot(1,4,4)
imagesc(Xe(1,:).*1000,Zvec.*1000,pIms(:,:,4),[10 100])
colormap('gray')
title('16 parallel foci')
xlabel('Xe (mm)')
ylabel('Depth (mm)')
axis image

%% Part 7 Apodization



