%  Make differene images for showing the effect of
%  focusing and apodization
%
%  This script assumes that the field_init procedure has been called
%
%  Example by Joergen Arendt Jensen, March 25, 1997.
%  constant F# added by Peter Munk, April 8, 1997
%  Version 1.2, August 13, 2007, JAJ: Printout changed
%
%  Generate the transducer apertures for send and receive

f0=3e6;                  %  Transducer center frequency [Hz]
fs=100e6;                %  Sampling frequency [Hz]
c=1540;                  %  Speed of sound [m/s]
lambda=c/f0;             %  Wavelength [m]
width=lambda;            %  Width of element
element_height=5/1000;   %  Height of element [m]
kerf=0.1/1000;           %  Kerf [m]
focus=[0 0 70]/1000;     %  Fixed focal point [m]
N_elements=128;          %  Number of physical elements
N_active=64;             %  Number of active elements
xmit_N_active=128;       %  Number of active transmit elements for constant F#
rec_N_active=128;	 %  Number of active receive elements for constant F#

%  Set the sampling frequency 

set_sampling(fs);

%  Generate aperture for emission

emit_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 1,focus);

%  Set the impulse response and excitation of the emit aperture

impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse (emit_aperture, impulse_response);

excitation=sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation (emit_aperture, excitation);

%  Generate aperture for reception

receive_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 1,focus);

%  Set the impulse response for the receive aperture

xdc_impulse (receive_aperture, impulse_response);

%   Load the computer phantom
[phantom_positions, phantom_amplitudes] = cyst_phantom(10000);

%   Do linear array imaging

no_lines=20;                         %  Number of lines in image
image_width=no_lines/1000;                 %  Size of image sector
d_x=image_width/no_lines;            %  Increment for image

%  Use 128 elements with 2 mm between receive focuses

mecr128a
%%

mk_img_tom
subplot(1,3,1)
imagesc(xx,zz(100:250),new_env(100:250,:),[-50 0])
axis image
colormap(gray)
title('Ones')

%%
mk_img_tom_forCNR

%% create a struct that saves means and stds for diff windows
% all the prep done in command line
win_results.uhamClut = uhamClut;
win_results.uhamCyst = uhamCyst;
win_results.stdhamClut = stdhamNCyst;
win_results.stdhamCyst = stdhamCyst;
win_results.uhanClut = uhanClut;
win_results.uhamCyst = uhanCyst;
win_results.stdhanClut = stdhanClut;
win_results.stdhanCyst = stdhanCyst;
win_results.uoneCyst = uonecyst;
win_results.uonenoise = uonenoise;
win_results.stdonecyst = stdonecyst;
win_results.stdonenoise = stdonenoise;

save('windowing_results.mat','win_results');




%%
fig = gcf;
new_env = fig.Children(3).Children.CData;
%%
subplot(1,3,3)
image(XX,ZZ,new_env)
colormap(gray(128))
axis('image')
axis([-10 10 10 120])
title('Hamming')

%% get signal and noise regions
cyst_mask = roipoly();
noise_mask = roipoly();




%% Part 8 Constant F#
% make a mask showing binary aperature for different depths
%calc number of elements to use
dim = size(image_data);
Fnumber_rec = 2;
rec_N_active_dyn=round(((1:n)/fn+min_sample/fs)*1540/2/(Fnumber_rec*(width+kerf)));
% construct 2d matrix


%%
%implement
f0=3e6;                  %  Transducer center frequency [Hz]
fs=100e6;                %  Sampling frequency [Hz]
c=1540;                  %  Speed of sound [m/s]
lambda=c/f0;             %  Wavelength [m]
width=lambda;            %  Width of element
element_height=5/1000;   %  Height of element [m]
kerf=0.1/1000;           %  Kerf [m]
focus=[0 0 70]/1000;     %  Fixed focal point [m]
N_elements=128;          %  Number of physical elements
N_active=64;             %  Number of active elements
xmit_N_active=128;       %  Number of active transmit elements for constant F#
rec_N_active=128;	 %  Number of active receive elements for constant F#

%  Set the sampling frequency 

set_sampling(fs);

%  Generate aperture for emission

emit_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 1,focus);

%  Set the impulse response and excitation of the emit aperture

impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse (emit_aperture, impulse_response);

excitation=sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation (emit_aperture, excitation);

%  Generate aperture for reception

receive_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 1,focus);

%  Set the impulse response for the receive aperture

xdc_impulse (receive_aperture, impulse_response);

%   Load the computer phantom

[phantom_positions, phantom_amplitudes] = pts_pha;

%   Do linear array imaging

no_lines=20;                         %  Number of lines in image
image_width=20/1000;                 %  Size of image sector
d_x=image_width/no_lines;            %  Increment for image



 %  Use 128 elements instead with apodization 
 % constant F# rec=2 and xmit=4

disp('Constant F# number focusing')
%  The different focal zones for transmit

xmit_focal_zones_center=[10 20 40 80]'/1000; % focuspoints
xmit_focal_zones=[0 15 30 60]'/1000; % start of zone [m]
Nft=max(size(xmit_focal_zones));
focus_times_xmit=xmit_focal_zones/1540;


%  Set a Hanning apodization on the receive aperture
%  Dynamic opening aperture is used.

Fnumber_xmit=4.0;
xmit_N_active_dyn=round(xmit_focal_zones_center./(Fnumber_xmit*(width+kerf)));

for ii=1:Nft
 if xmit_N_active_dyn(ii)>xmit_N_active 
   xmit_N_active_dyn(ii)=xmit_N_active; end
 xmit_N_pre_dyn(ii) = ceil(xmit_N_active/2  - xmit_N_active_dyn(ii)/2);
 xmit_N_post_dyn(ii) = xmit_N_active - xmit_N_pre_dyn(ii) - xmit_N_active_dyn(ii);
 xmit_apo=(hanning(xmit_N_active_dyn(ii)))';
 xmit_apo_matrix_sub(ii,:)=[zeros(1,xmit_N_pre_dyn(ii)) xmit_apo zeros(1,xmit_N_post_dyn(ii))];
end
apo_zero=(N_elements-xmit_N_active)/2;
xmit_apo_matrix=[zeros(Nft,apo_zero) xmit_apo_matrix_sub zeros(Nft,apo_zero)];

%**********************************************************
%                      RECEIVE
%**********************************************************

%  The different focal zones for receive

%  Set the different focal zones for reception

rec_zone_start=30/1000;
rec_zone_stop=150/1000;
rec_zone_size=20/1000;

rec_focal_zones_center=[rec_zone_start:rec_zone_size:rec_zone_stop]';
rec_focal_zones=rec_focal_zones_center-0.5*rec_zone_size;
Nfr=max(size(rec_focal_zones));
focus_times_rec=rec_focal_zones/1540;


%  Set a Hanning apodization on the receive aperture
%  Dynamic opening aperture is used.

Fnumber_rec=2.0;
rec_N_active_dyn=round(rec_focal_zones_center./(Fnumber_rec*(width+kerf)));
%%
for ii=1:Nfr
 if rec_N_active_dyn(ii)>rec_N_active 
   rec_N_active_dyn(ii)=rec_N_active; end
 rec_N_pre_dyn(ii) = ceil(rec_N_active/2  - rec_N_active_dyn(ii)/2);
 rec_N_post_dyn(ii) = rec_N_active - rec_N_pre_dyn(ii) - rec_N_active_dyn(ii);
 rec_apo=(hanning(rec_N_active_dyn(ii)))';
 rec_apo_matrix_sub(ii,:)=[zeros(1,rec_N_pre_dyn(ii)) rec_apo zeros(1,rec_N_post_dyn(ii))];
end
apo_zero=(N_elements-rec_N_active)/2;
rec_apo_matrix=[zeros(Nfr,apo_zero) rec_apo_matrix_sub zeros(Nfr,apo_zero)];

%**********************************************************
%                      MAKE IMAGE
%**********************************************************

x= -image_width/2;
image_data=0;
for i=1:no_lines
 % Set the focus for this direction (aperture times points(x,y,z))

  xdc_center_focus (emit_aperture,[x 0 0]);
  xdc_focus (emit_aperture, focus_times_xmit, [x*ones(Nft,1), zeros(Nft,1), xmit_focal_zones_center]);
  xdc_center_focus (receive_aperture,[x 0 0]);
  xdc_focus (receive_aperture, focus_times_rec, [x*ones(Nfr,1), zeros(Nfr,1), rec_focal_zones_center]);

   %  Calculate the apodization 

  xdc_apodization (emit_aperture, focus_times_xmit, xmit_apo_matrix);
  xdc_apodization (receive_aperture, focus_times_rec , rec_apo_matrix);

  %   Calculate the received response

  [v, t1]=calc_scat(emit_aperture, receive_aperture, phantom_positions, phantom_amplitudes);

  %  Store the result

  image_data(1:max(size(v)),i)=v;
  times(i) = t1;

  %  Steer in another angle

  x = x + d_x;
end
mk_img
title('Constant F#')
%%
% show binary image of aperature mask
imagesc(1:128,rec_focal_zones_center.*1000, rec_apo_matrix>.07)
xlabel('Elements')
ylabel('Focal Depth (mm) (I used 7 discrete depths)')
title('F# mask')

%% apodized, dynamic focus, but not constant f#

figure
mecr128a
mk_img
title('Hamming')

