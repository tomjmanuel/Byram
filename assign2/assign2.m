%% Assignment 2 Tom Manuel
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


%%





