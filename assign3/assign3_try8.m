%% Synthetic aperature
% try6 has mask 
% adapted from parallel imaging code...
% with parallel imaging
% each foci will have a slightly different delay profile
% so calculate a Trx for each
close all
clear all
load('pointTargetData.mat')

pff = 256; % parallel focus factor (2,4,8,16)
data = veraStrct.data;


%First I will create a grid that represents the x and z coordinates of each
%pixel in physical space
Nx = 2.*veraStrct.numElementsPerXmt;
dx = 1E-3*veraStrct.XMTspacingMM;

%get Nx,dx,Nz, and dz
foo = size(data);
Nz = foo(1);
t0 = veraStrct.timeZero -1; % nPts to throw away 
Nz = foo(1)-t0;
data = data(t0+1:end,:,:);
%foo = size(data);
foo = [Nz Nx Nx];

% dz = [m/s]*[s/sample] = c * fs^-1 
c = 1540; %m/s
dz = .5 * c / 20E6;

% build coordinate matrices (centered relative to element positions)
xx = linspace(-(Nx-1)*dx/2,(Nx-1)*dx/2,Nx);
Xe = repmat(xx,[Nz 1]);
Ze = repmat(dz.*linspace(0,Nz-1,Nz)',[1,Nx]);
t = 2.*Ze(:,1)./c;

% create lateral foci vector (aline position)
Xf = xx;
Zf = .04;

%% Calculate Ttx for all points in grid, relative to the transmit beam location

% Ttx is (zf - radius_p0)/c;
% Calculate radius P0; (distance between each point and the focus
%radMat = sqrt((transmitBeamLoc-Xe).^2 + (Zf-Ze).^2);
radMat = sqrt((Xe).^2 + (Zf-Ze).^2);

%zflag is -1 before focus, positive one after
zflag = ones(size(Ze));
zflag(Ze<Zf)=-1;
Ttx = (Zf + zflag.*radMat)/c;

%% Calculate Trx for each aline that goes with our one beam
% this delay corresponds to distance from 
Trx = zeros(foo);
for i=1:pff
    Trx(:,:,i)=(1/c).*sqrt((Xf(i)-Xe).^2 + Ze.^2);
end

%% Combine Trx and Ttx
%[Z channel aline]
Ttot = zeros(foo);
for i=1:pff
    Ttot(:,:,i) = repmat(Ttx(:,i),[1 Nx]) + Trx(:,:,i);
end

%% make intensity mask that blocks out things not in beam path
%estimate with a parabola
bw = dx.*3; %beam width at focus
p2 = (((Nx/2+1).*dx)/2-bw)/Zf^2.*(Ze-Zf).^2 + bw;
mask = abs(Xe)<p2;

%% apply delays
ddata = zeros(foo);
%h = waitbar(0,'computing');
finalIm = zeros([Nz Nx]);
%for bb = 1:128
for bb=50:80
    beam = data(:,:,bb);
    beam = padarray(beam,[0 64]);
    %
    for i=1:pff %loop through alines
        for j=1:Nx % loop through channel data
            ddata(:,j,i) = interp1(t,beam(:,j),Ttot(:,j,i));
        end
    end
    
    imf = mask.*squeeze(sum(ddata,2));
    
    %imagesc(imf,[10 100])
    %pause(.1)
    
    finalIm(:,bb:bb+128-1)=finalIm(:,bb:bb+128-1)+imf(:,65:128+64);
    %waitbar(bb/128,h);
    bb
end
%close(h);


%%


imf = finalIm;
imf(isnan(imf))=-1000; %interp1 drops some Nans in here
imf = 20.*log10(abs(hilbert(imf)));
%%
% image of delayed and summed data
Zvec = dz.*linspace(0,Nz-1,Nz); % z axis for plotting time
figure
imagesc(xx.*1000,Zvec.*1000,imf,[15 120])
colormap('gray')
title('Parallel Imaging')
xlabel('Xe (mm)')
ylabel('Depth (mm)')
axis image
