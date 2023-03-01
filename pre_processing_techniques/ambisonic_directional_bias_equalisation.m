function [SH_ambisonic_binaural_decoder_DBE,VL_hrirs_DBE] = ambisonic_directional_bias_equalisation(ambisonic_order,SH_ambisonic_binaural_decoder,bias_factor,bias_direction,VL_hrirs,Fs,dbe_plot_flag)
%{
Function to augment the frequency spectra of virtual loudspeaker HRIRs for
Ambisonic reproduction, in order to improve the Ambisonic timbral accuracy
in the specified direction of bias.

Thomas McKenzie, University of York, 2019.

References:
"Directional bias equalisation of first order binaural Ambisonic rendering,"
in AES Conference on Audio for Virtual and Augmented Reality, Redmond,
USA, 2018, T. McKenzie, D. T. Murphy, and G. Kearney.

"Towards a perceptually optimal bias factor for directional bias equalisa-
tion of binaural Ambisonic rendering," in EAA Spatial Audio Signal Processing
Symposium, (pp. 97{102), Paris, France, 2019, T. McKenzie, D. T. Murphy, and G.
Kearney.
%}

%%
bias_direction_azimuth = bias_direction(1);
bias_direction_elevation = bias_direction(2);

window_size = length(SH_ambisonic_binaural_decoder(1,:,1));

%% Get points on a sphere and encode them to Ambisonics

points_in_RMS_calculation = 1000;
[azVec, elVec] = fibonacci2azi_ele(points_in_RMS_calculation);
[azVec, elVec] = bias_quadrature(azVec,elVec,bias_direction_azimuth...
    ,bias_direction_elevation,bias_factor,dbe_plot_flag);

azi_ele_matrix = [azVec, elVec];
ambisonic_renders = zeros(points_in_RMS_calculation,window_size,2);

for i = 1:length(azVec)
    encodedAmbisonicGains = encodeHOA_N3D(ambisonic_order,1,azi_ele_matrix(i,:));
    ambisonicHRIRs(1,:) = encodedAmbisonicGains * SH_ambisonic_binaural_decoder(:,:,1);
    ambisonicHRIRs(2,:) = encodedAmbisonicGains * SH_ambisonic_binaural_decoder(:,:,2);
    ambisonic_renders(i,:,:) = ambisonicHRIRs';
end

% Diffuse Field Equalisation parameters
octave_band_smoothing_fraction = 8; % Octave band smoothing (0 = off, 1 = Octave, 2 = 1/2 Octave etc) 4 = 1/4 octave smoothing.
inversion_range = [2 20000]; % Range for inversion
regularisation_parameters = [30 20]; % In band and out of band regularisation parameters (dB)

hrir_renders = zeros(points_in_RMS_calculation, (window_size-1), 2); % Initialise output HRIRs
L_AVG = zeros(window_size,1); % Initialise average left ear response
R_AVG = zeros(window_size,1); % Initialise average right ear response

for i = 1:points_in_RMS_calculation
    hrir_renders(i,:,1) = ambisonic_renders(i,1:(window_size-1),1);
    hrir_renders(i,:,2) = ambisonic_renders(i,1:(window_size-1),2);
end

% Contribution of HRIR to average response is NOT based on solid angle -
% every direction weighted the same
for i = 1:points_in_RMS_calculation
    HRIR_L = fft(hrir_renders(i,:,1)', window_size);
    HRIR_R = fft(hrir_renders(i,:,2)', window_size);
    
    % Average response *NOT* based on the solid angle:
    L_AVG = L_AVG + abs(HRIR_L).^2;
    R_AVG = R_AVG + abs(HRIR_R).^2;
end

L_AVG = L_AVG / points_in_RMS_calculation;
R_AVG = R_AVG / points_in_RMS_calculation;
L_AVG = sqrt(L_AVG); R_AVG = sqrt(R_AVG);

dbq_avg_L = rotate_vect(real(ifft(L_AVG)),window_size/2);
dbq_avg_R = rotate_vect(real(ifft(R_AVG)),window_size/2);

% Compute Inverse Filters
[ih_L]=invFIR('linphase',dbq_avg_L,window_size,octave_band_smoothing_fraction...
    ,window_size,inversion_range,regularisation_parameters,1, Fs);
[ih_R]=invFIR('linphase',dbq_avg_R,window_size,octave_band_smoothing_fraction...
    ,window_size,inversion_range,regularisation_parameters,1, Fs);

%% Directionally Biased Quadrature frontal HRIR weighting

bias_factor_scaled = 1 - exp(-0.1*(bias_factor-1));
biasDirHRIR = audioread(char(strcat('hrirs/azi_',num2str(bias_direction_azimuth)...
    ,'_ele_',num2str(bias_direction_elevation),'_DFC.wav')));
biasDirHRIR_compressed = compress_IR(biasDirHRIR,window_size,bias_factor_scaled);

ih_L2 = conv(biasDirHRIR_compressed(:,1),ih_L);
ih_R2 = conv(biasDirHRIR_compressed(:,2),ih_R);
ih_L = ih_L2(round(length(ih_L2)*3/8):(length(ih_L2)-round(length(ih_L2)*3/8)));
ih_R = ih_R2(round(length(ih_R2)*3/8):(length(ih_R2)-round(length(ih_R2)*3/8)));

%% Convolve original loudspeaker HRIRs with inverse filters, truncate and window

fade_windows = { @(N)(hanning(N).^2) @(N)(hanning(N).^2) };
ih_L = fade_samples(ih_L,[50 50], fade_windows);
ih_R = fade_samples(ih_R,[50 50], fade_windows);

% initialise new decoder
SH_ambisonic_binaural_decoder_DBE = zeros(length(SH_ambisonic_binaural_decoder(:,1,1))...
    ,length(SH_ambisonic_binaural_decoder(1,:,1)),length(SH_ambisonic_binaural_decoder(1,1,:)));

for i = 1:length(SH_ambisonic_binaural_decoder(:,1,1))
    SH_ambisonic_binaural_decoderEQ1(i,:,1) = conv(ih_L,SH_ambisonic_binaural_decoder(i,:,1));
    SH_ambisonic_binaural_decoderEQ1(i,:,2) = conv(ih_R,SH_ambisonic_binaural_decoder(i,:,2));

    SH_ambisonic_binaural_decoder_DBE(i,:,:) = SH_ambisonic_binaural_decoderEQ1(i,(length(ih_L)/2)+1:(length(SH_ambisonic_binaural_decoderEQ1(1,:,1)) - (length(ih_L)/2))+1,:); % THE +1 ensures a delay is not added of one sample (ttm 13.6.2018)
end

VL_hrirs_DBE = zeros(length(VL_hrirs(:,1,1))...
    ,length(VL_hrirs(1,:,1)),length(VL_hrirs(1,1,:)));
for i = 1:length(VL_hrirs(1,1,:))
    VL_hrirs_DBE1(1,:,i) = conv(ih_L,VL_hrirs(1,:,i));
    VL_hrirs_DBE1(2,:,i) = conv(ih_R,VL_hrirs(2,:,i));
    
    VL_hrirs_DBE(:,:,i) = VL_hrirs_DBE1(:...
        ,(length(ih_L)/2)+1:(length(VL_hrirs_DBE1(1,:,1)) - ...
        (length(ih_L)/2))+1,i); % THE +1 ensures a delay is not added of one sample (ttm 13.6.2018)
end

%% Plot results?

if dbe_plot_flag == 1
    comp_L = conv(ih_L, dbq_avg_L);
    plotColour = get(gca,'colororder');
    figure
    colour_frequency_plot(dbq_avg_L,Fs, window_size*4,  '-',plotColour(1,:), 1.8,[70 20000]);
    hold on
    colour_frequency_plot(comp_L,Fs, window_size*4, '-.',plotColour(2,:), 1.8,[70 20000]);
    colour_frequency_plot(ih_L,Fs, window_size*4, '-',[0 0 0], 1.4,[70 20000]);
    legend({'DBQ RMS Response' ,...%(L)', 'Diffuse Field Response (R)', ...
        'Directional Bias HRTF',...% (L)', 'Inverse Filter (R)',   ...
        'Equalisation Filter'}, 'location', 'southwest','FontSize',14);
    
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    set(gcf, 'Color', 'w');
    set(gca, 'fontsize', 14);
    xlim([70 20000]);
    ylim([-25 20]);
end
disp(strcat('Ambisonic DBE complete! Bias factor = ',num2str(bias_factor)));
end

function [azimuthOut, elevationOut] = bias_quadrature(azimuthIn,elevationIn,biasDirAz,biasDirEl,biasFactor,plotFlag)
% Input arrays of azimuth and elevation angles and bias the quadrature
% points towards a specified direction (biasDirAz and biasDirEl). If
% biasFactor is 1 then the output will be the same as the input. Returns
% arrays of azimuth and elevation in same format as input.
%
% ========================
% Inputs:
%     azimuthIn:        Vector of azimuth angles (0 to 360 degrees)
%
%     elevationIn:      Vector of elevation angles (-90 to 90 degrees)
%
%     biasDirAz:        Azimuth angle for bias direction (degrees)
%
%     biasDirEl:        Elevation angle for bias direction (degrees)
%
%     biasFactor:       Enter number > 1. Higher number --> higher biasing
%
%     plotFlag:         Plot voronoi sphere of output
%
%
% ========================
% Outputs:
%     azimuthOut:       Vector of azimuth angles (0 to 360 degrees)
%
%     elevationOut:     Vector of elevation angles (-90 to 90 degrees)
%
% ///////////// Tom McKenzie, University of York, 2018 \\\\\\\\\\\\\\\
%
% ========================
%

% convert spherical coordinates to cartesian coordinates
[x,y,z] = sph2cart(azimuthIn*pi/180, elevationIn*pi/180, ones(size(azimuthIn)));

% bias upwards (toward the north pole)
z1 = z + 1;
z2 = z1 * biasFactor;
z3 = z2 - 1;

% re-convert new cartesian coordinates to spherical coordinates
[azVec3, elVec3,~] = cart2sph( x,y,z3);
azVec4 = rad2deg(azVec3);

elevationOut1 = (rad2deg(elVec3))';
azimuthOut1 = zeros(1,length(azimuthIn));
for i = 1:length(azVec4)
    if azVec4(i) < 0
        azimuthOut1(i) = azVec4(i) + 360;
    else
        azimuthOut1(i) = azVec4(i);
    end
end

% rotate the entire quadrature (now biased) to the desired bias direction
azimuthOut = zeros(length(azimuthIn),0);
elevationOut = zeros(length(elevationIn),0);
for i = 1:length(azimuthIn)
    [mediumAzEl] = rotate_azi_ele_alt([azimuthOut1(i) elevationOut1(i)],[0 -90]);
    [newAzEl] =  rotate_azi_ele_alt(mediumAzEl,[biasDirAz biasDirEl]);
    
    azimuthOut(i) = newAzEl(1);
    elevationOut(i) = newAzEl(2);
end

% plot if requested
if plotFlag == 1
    GetVoronoiPlotandSolidAng(azimuthOut, elevationOut, 1);
    colormap(flipud(winter))
end

% rotate vectors
azimuthOut = azimuthOut';
elevationOut = elevationOut';
end

function [ azi,ele ] = fibonacci2azi_ele( noOfPoints )
% fibonacci2azi_ele()
%
% Converts fibonacci spiral quadrature to a matrix of azimuth and elevation.
% Azimuth angles 0 to 360 (going counterclockwise)
% Elevation angles -90 to 90 (going upwards)

[xyz] = sphere_fibonacci_grid_points(noOfPoints);

% initiate output arrays
azi = zeros(length(xyz),1); ele = zeros(length(xyz),1);

for i = 1:length(xyz)
    [azi(i),ele(i)] = xyz_to_tp_With_P_origin_at_0(xyz(i,1),xyz(i,2),xyz(i,3));
end
end

function xg = sphere_fibonacci_grid_points ( ng )

%*****************************************************************************80
%
%% SPHERE_FIBONACCI_GRID_POINTS: Fibonacci spiral gridpoints on a sphere.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    27 April 2015
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Richard Swinbank, James Purser,
%    Fibonacci grids: A novel approach to global modelling,
%    Quarterly Journal of the Royal Meteorological Society,
%    Volume 132, Number 619, July 2006 Part B, pages 1769-1793.
%
%  Parameters:
%
%    Input, integer NG, the number of points.
%
%    Output, real XG(N,3), the grid points.
%
phi = ( 1.0 + sqrt ( 5.0 ) ) / 2.0;

i = ( - ( ng - 1 ) : 2 : ( ng - 1 ) )';
theta = 2 * pi * i / phi;
sphi = i / ng;
cphi = sqrt ( ( ng + i ) .* ( ng - i ) ) / ng;

xg = zeros ( ng, 3 );

xg(1:ng,1) = cphi .* sin ( theta );
xg(1:ng,2) = cphi .* cos ( theta );
xg(1:ng,3) = sphi;

return
end

function [compressedIR] = compress_IR(inputIR,Nfft,biasFactor)

% ------------------------------------------------------------------------
% function - compress_IR()
% Enter a biasFactor value less than 1 to compress the HRIR. a value of 0
% will mean a flat spectrum at the mean amplitude of the original HRIR and
% a vaule of 1 will mean the same as the original HRIR.
% ---------------
%
% inputImpulseResponse  - mono or stereo impulse response (column vector)
%
% Nfft                  - FFT length for calculating compressed FIR
%
% biasFactor            - number between 0 and 1 (0 flattens the response
%                         completely, 1 returns the exact original IR)
%
%
% Thomas McKenzie, University of York, 20 / 4 / 2018

% calculate compressed filter
H = abs(fft(inputIR(:,1),Nfft));
H2 = (H*biasFactor)+((1-biasFactor)*(mean(H)));

compressedIR=circshift(ifft(H2,'symmetric'),Nfft/2);
compressedIR=compressedIR(end/2-Nfft/2+1:end/2+Nfft/2); % truncation to length L

% 2-channel case
%-------------------------------------------------
if size(inputIR,2)==2
    H = abs(fft(inputIR(:,2),Nfft));
    H2 = (H*biasFactor)+((1-biasFactor)*(mean(H)));
    
    compressedIRr=circshift(ifft(H2,'symmetric'),Nfft/2);
    compressedIRr=compressedIRr(end/2-Nfft/2+1:end/2+Nfft/2);
    compressedIR=[compressedIR compressedIRr];
end
end

function new_azi_ele = rotate_azi_ele_alt(old_azi_ele,rotateby_azi_ele)
% Alteration of rotate() (by Matlab) by Thomas McKenzie, University of York
% ROTATE() -- Copyright 1984-2017 The MathWorks, Inc.

% first rotate by azimuth (easy)
old_azi_ele(1) = old_azi_ele(1) + rotateby_azi_ele(1);
if old_azi_ele(1) >= 360
    old_azi_ele(1) = old_azi_ele(1) - 360;
end
if old_azi_ele(1) < 0
    old_azi_ele(1) = old_azi_ele(1) + 360;
end

[oldXYZ(1),oldXYZ(2),oldXYZ(3)] = sph2cart(deg2rad(old_azi_ele(1)),deg2rad(old_azi_ele(2)),1);

% find unit vector for axis of rotation
theta = pi*(rotateby_azi_ele(1)-90)/180;
phi = pi*0/180;
u = [cos(phi)*cos(theta); cos(phi)*sin(theta); sin(phi)];

alph = (rotateby_azi_ele(2))*pi/180;
cos_a = cos(alph);
sin_a = sin(alph);
ver_a = 1 - cos_a;
x = u(1);
y = u(2);
z = u(3);
rot = [cos_a+x^2*ver_a x*y*ver_a-z*sin_a x*z*ver_a+y*sin_a; ...
    x*y*ver_a+z*sin_a cos_a+y^2*ver_a y*z*ver_a-x*sin_a; ...
    x*z*ver_a-y*sin_a y*z*ver_a+x*sin_a cos_a+z^2*ver_a]';

newXYZ = oldXYZ*rot;

[newAzEl1(1),newAzEl1(2),~] = cart2sph(newXYZ(1),newXYZ(2),newXYZ(3));
new_azi_ele(1) = rad2deg(newAzEl1(1));
new_azi_ele(2) = rad2deg(newAzEl1(2));
end


