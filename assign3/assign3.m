%% Implement Synthetic aperature as described in 
% A Study of Synthetic-Aperture Imaging with
% Virtual Source Elements in B-Mode
% Ultrasound Imaging Systems
% Tom Manuel, 10/11/19 Assign 3

clear all
close all
load('pointTargetData.mat')

data = veraStrct.data;
t0 = veraStrct.timeZero -1; % nPts to throw away 
data = data(t0+1:end,:,:);
% create a grid that represents the x and z coordinates of each
% pixel in physical space
Nx =  2*veraStrct.numElementsPerXmt;
dx = 1E-3*veraStrct.XMTspacingMM;
fs = 1E6*veraStrct.samplingRateMHz;
%get Nx,dx,Nz, and dz
foo = size(data);
Nz = foo(1);

Nz = foo(1)-t0;
data = data(t0+1:end,:,:);
foo = size(data);
% dz = [m/s]*[s/sample] = c * fs^-1 
c = 1540; %m/s
v=c;
dz = .5 * c / fs;

% build coordinate matrices
Xef = repmat(dx.*linspace(-Nx/2,Nx/2,Nx),[Nz 1]); %x coords
Zef = repmat(dz.*linspace(0,Nz-1,Nz)',[1,Nx]); %y coords

%% implement beamforming for single beam first
% calculate delay profile for each point in the image
frame=64;
bdata = data(:,:,frame); % one beam
bdata = padarray(bdata,[0 frame]);

% calc delays
tdmat = zeros(Nx,Nx,Nz); %[alines, channel, time]
% get transmit delay for desired channel
t = 2.*Zef(:,1)./c; %round trip time between center element and zp;
zf = 4E-2; %transmit focus at 4cm
%xf is distance from desired a line to center of beam
bc = (frame+64)*dx; %beam center
al = i*dx; % aline location

Xe1d= Xef(1,:)';

%flag to say if above or below focus
zflag = ones(Nz,1);
zflag(1:round(zf*dz/2))=-1;

%synthesized a lines
for i=1:256
    al = i*dx; % aline location
    xf = abs(bc -al);
    % trasnd holds transmit delay for an a line off centered from beam 
    % eqn [1]
    transd = t/2 - (zf + zflag.*sqrt(xf^2 + (zf - v.*t/2).^2))/v;
    
    % now get eqn [5] (receive delay)
    %channels
    for j=1:256
        %xr is distance between a line and element
        xr = abs(Xe1d(j)-al);
        foo1 = sqrt((zf - (v.*t/2)).^2+xf^2);
        foo2 = sqrt((v.*t/2).^2+xr^2);
        recd = t - ((zf-foo1)/v + foo2/v); % receive delay (s)
        tdmat(i,j,:) = transd + recd;
        
    end
end

%cast tdmat from s to samp
tdmat = tdmat.*fs;
tdmat = tdmat - min(tdmat(:));


%%
%create time vector
t = linspace(0,Nz,Nz)';

% delay data
beams = zeros(Nx,Nx,Nz);

%for each a line
for i=1:Nx
    %for each xr
    for j=1:Nx
        beams(i,j,:) = interp1(t,bdata(:,j),squeeze(tdmat(i,j,:)));
    end
end

%%
imf = squeeze(sum(beams,2))';

imf(isnan(imf))=-1000; %interp1 drops some Nans in here
% compress
imf = 20.*log10(abs(hilbert(imf)));

% image of delayed and summed data
Zvec = dz.*linspace(0,Nz-1,Nz); % z axis for plotting time
figure
imagesc(Xe1d.*1000,Zvec.*1000,imf,[10 100])
colormap('gray')
title('Parallel Imaging')
xlabel('Xe (mm)')
ylabel('Depth (mm)')
axis image

%%









