%% DEFINE ARRAY L12
c = 1540;                   % Speed of sound
f0 = 7.813e6;                 % Transducer center frequency [Hz]
fs = 31.25*1e6;                 % Sampling frequency [Hz]
% fs = 250e6;                 % Sampling frequency [Hz]
lambda = c/f0;              % Wavelength
element_height = 5/1000;   % Height of element [m] (elevation direction)
pitch = 0.2/1000;         % Distance between element centers
element_width = 0.17/1000; % Element width [m] (azimuth direction)
kerf = pitch-element_width;          % Width of fill material between the ceramic elements
Rfocus = 20/1000;           % Elevation lens focus (or radius of curvature, ROC) 
focus = [0 0 inf]/1000;      % Fixed emitter focal point [m] (irrelevant for single element transducer)
N_elements = 64;            % Number of physical elements in array
N_sub_x = 10;                % Element sub division in x-direction
N_sub_y = 10;               % Element subdivision in y-direction
apodTx = 0;                 % Transmit apodization. 0:boxcar, 1:Hanning, 2:Cosine-tapered 0.3
apodRx = 1;                 % Receive apodization. 0:boxcar, 1:Hanning, 2:Cosine-tapered 0.3
dynamic_receive_focus = 0;  % Enable dynamic receive focusing.
simType = 'txrx';
% Optional measured-waveform assets and external dependencies can be kept
% outside the repository. Point VRUSI_ASSETS_ROOT to that local directory.
assets_root = getenv('VRUSI_ASSETS_ROOT');
if ~isempty(assets_root)
    addpath(genpath(assets_root));
end

param = getparam('L12-3v');
param.fc = f0;
param.kerf = kerf;
param.width = element_width;
param.pitch = pitch;
param.Nelements = N_elements;
param.height = element_height;
param.focus = Rfocus;
param.fs = fs; % sampling frequency
param.TXapodization = [ones(1,N_elements)];

%% SET THE IMPULSE RESPONSE AND EXCITATION OF THE TRANSMIT AND RECEIVE APERTURE 
% t_ir = -2/f0:1/(fs):2/f0;
% Bw = 0.8;
% impulse_response=gauspuls(t_ir,f0,Bw);

load('tv.mat');
impulse_response = downsample(tv1(:,2),16*8);
impulse_response = impulse_response(35:55,1);
% impulse_response = impulse_response(70:150,1);
impulse_response = rescale(impulse_response,-1,1)';
impulse_response = impulse_response - mean(impulse_response);
set_sampling(fs);

%% 
% ex_periods = 2.0;
% t_ex=(0:1/(fs):ex_periods/f0);
% excitation=square(2*pi*f0*t_ex);
load('TW_L12_067DC_2cycles.mat');
excitation = downsample(TW.Wvfm1Wy',8);

%% DEFINE MEASUREMENT POINTS
wavelengthh = lambda;%in m
pdelta = [0.519480519480520,0,0.500000000000000];% in wavelength
porigin = [-32.727272727272734,0,10.146753246753248];% in wavelength
measDepthStart = porigin(1,3)*wavelengthh;    % Start depth along z-axis to place measurement points in mm
measDepthEnd = 1.522012987012987e+02*wavelengthh;    % End depth along z-axis to place measurement points
xStart = porigin(1,1)*wavelengthh;      % Start position of measurement points in x direction, in mm
xEnd = -xStart;         % End position of measurement points in x direction

Nmpx = N_elements*2;
Nmpz = 285;

mx = linspace(xStart,xEnd,Nmpx)';
my = zeros(Nmpx*Nmpz, 1);
mz = linspace(measDepthStart,measDepthEnd,Nmpz)';

[X,Z] = meshgrid(mx,mz);
measurement_points = [X(:),my,Z(:)];
