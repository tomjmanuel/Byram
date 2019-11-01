
%% Prob1 assign 3
% Tom Manuel
% Load in the data

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

%% Filter Data 
% select parameters
fs = 1E6.*veraStrct.samplingRateMHz; %(sam/s)

Nf = fs/2; %nyquist freq (sam/s)

Fc = 1E6.*veraStrct.frequencyMHz; %center freq (cyc/s)
Fcs = Fc/fs; % center freq (cyc/sam)

%
aline = data(:,64,64);
close all
bandwidth = .08;
%cf = .55; %center frequency
cf = Fcs;
lf = cf -bandwidth;
hf = cf+bandwidth;

% create filter ([order, window]) with window being relative to samp rate
b = fir1(20,[lf hf]);
dfilt = filter(b,1,data);

%%


% woah do you have to guess c to get dz??
% dz = [m/s]*[s/sample] = c * fs^-1 
c = 1540; %m/s
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

