%  Compress the data to show 60 dB of
%  dynamic range
%
%  Example by Joergen Arendt Jensen, Nov. 28, 1995.

%  Adjust the data in time and display it as
%  a gray scale image
%
%  Version 1.1, April 1, 1998, Joergen Arendt Jensen
%  Version 1.2, August 13, 2007, JAJ: 
%       Printout changed
%       Display to a 60 dB dynamic range



D=1;             %  Sampling frequency decimation factor
% the decimation factor downsamples in time domain

% first setup 'n' to add appropriate samples to beginning of image
min_sample=9/1000/c;
max_sample=max([max(times)*fs 2*0.121/c*fs]);
[n,m]=size(image_data);

n=n+(max_sample-min_sample);
env=zeros(round(n/D),no_lines);

for i=1:no_lines
  % very sneak way of adding zeros to beginning of image_data one line at a
  % time
  % [ A; B] concatenates arrays
  rf_env=abs(hilbert([zeros(round(times(i)*fs-min_sample),1); image_data(:,i)]));
 % Here downsampling is performed
  rf_env=rf_env(1:D:max(size(rf_env)));
  
  % store column in envelope image (2d array)
  env(1:max(size(rf_env)),i)=rf_env;

end

%  Do logarithmic compression to a 50 dB dynamic range

log_env=20*log10(env);
log_env=log_env-max(max(log_env));

% set lower cutoff at 200 dB
log_env(log_env<-200)=-200;

%  Make an interpolated image
ID=1;
[n,m]=size(log_env);
new_env=zeros(n,m*ID);
for i=1:n
  new_env(i,:)=interp(log_env(i,:),ID);
  end
[n,m]=size(new_env);
  
% now plot with appropriate axis
% we need n points that range from 0 to Zmm and m points for 0 to Xmm
% what is Zmm? dpdentent on fs, speed of sound, minimum sample from times,
fn=fs/D;
zz =((1:n)/fn+min_sample/fs)*1540/2*1000;

% Xmm is dependent on no_lines, dx, and interpolation factor
xx = ((1:(ID*no_lines))*d_x/ID-no_lines*d_x/2)*1000;

figure
imagesc(xx,zz,new_env,[-50 0])
%axis image
colormap(gray)
