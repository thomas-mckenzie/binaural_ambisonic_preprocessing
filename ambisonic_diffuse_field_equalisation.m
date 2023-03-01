function [SH_ambisonic_binaural_decoder_DFE,VL_hrirs_DFE,df_avg_L,df_avg_R] = ambisonic_diffuse_field_equalisation(ambisonic_order,SH_ambisonic_binaural_decoder,VL_hrirs,Fs,dfe_plot_flag)
%{
Function to augment the frequency spectra of virtual loudspeaker HRIRs for
Ambisonic reproduction, in order to improve the Ambisonic timbral accuracy
in the diffuse-field.

Thomas McKenzie, University of York, 2019.

References:
"Diffuse-field equalisation of binaural Ambisonic rendering," Applied Sci-
ences, vol. 8, no. 10, 2018, T. McKenzie, D. T. Murphy, and G. Kearney.
%}

%%
window_size = length(SH_ambisonic_binaural_decoder(1,:,1));

%% Get points on a sphere and encode them to Ambisonics

pointsOnSphere = 240;
[azVec, elVec] = tdesign2azi_ele( pointsOnSphere );
azi_ele_matrix = [azVec, elVec];
azVec = azVec'; elVec = elVec';

ambisonic_renders = zeros(pointsOnSphere,window_size,2);
for i = 1:length(azVec)
    encodedAmbisonicGains = encodeHOA_N3D(ambisonic_order,1,azi_ele_matrix(i,:));
    ambisonicHRIR(1,:) = encodedAmbisonicGains * SH_ambisonic_binaural_decoder(:,:,1);
    ambisonicHRIR(2,:) = encodedAmbisonicGains * SH_ambisonic_binaural_decoder(:,:,2);
    ambisonic_renders(i,:,:) = ambisonicHRIR';
end

% Diffuse Field Equalisation parameters
octave_band_smoothing_fraction = 4; % Octave band smoothing (0 = off, 1 = Octave, 2 = 1/2 Octave etc) 4 = 1/4 octave smoothing.
inversion_range = [2 20000]; % Range for inversion
regularisation_parameters = [30 20]; % In band and out of band regularisation parameters (dB)

hrir_renders = zeros(pointsOnSphere, (window_size-1), 2); % Initialise output HRIRs
for i = 1:pointsOnSphere
    hrir_renders(i,:,1) = ambisonic_renders(i,1:(window_size-1),1);
    hrir_renders(i,:,2) = ambisonic_renders(i,1:(window_size-1),2);
end

% Get weights based on solid angle
solid_angle_weighting = GetVoronoiPlotandSolidAng(azVec, elVec, 0); % change last variable to 1 to print out voronoi sphere plot
solid_angle_weighting = solid_angle_weighting/(sum(solid_angle_weighting));

L_AVG = zeros(window_size,1); % Initialise average left ear response
R_AVG = zeros(window_size,1); % Initialise average right ear response
% Contribution of HRIR to average response is dependent on solid angle
for i = 1:pointsOnSphere
    HRIR_L = fft(hrir_renders(i,:,1)', window_size);
    HRIR_R = fft(hrir_renders(i,:,2)', window_size);
    
    % Average response based on the solid angle:
    L_AVG = L_AVG + solid_angle_weighting(i)*abs(HRIR_L).^2;
    R_AVG = R_AVG + solid_angle_weighting(i)*abs(HRIR_R).^2;
end

L_AVG = sqrt(L_AVG); R_AVG = sqrt(R_AVG); % square root in freq domain
df_avg_L = rotate_vect(real(ifft(L_AVG)),window_size/2); % convert to time domain
df_avg_R = rotate_vect(real(ifft(R_AVG)),window_size/2);

% Compute Inverse Filters
[ih_L]=invFIR('linphase',df_avg_L,window_size,octave_band_smoothing_fraction...
    ,window_size,inversion_range,regularisation_parameters,1, Fs);
[ih_R]=invFIR('linphase',df_avg_R,window_size,octave_band_smoothing_fraction...
    ,window_size,inversion_range,regularisation_parameters,1, Fs);

%% Convolve original loudspeaker HRIRs with inverse filters, truncate and window

fade_windows = { @(N)(hanning(N).^2) @(N)(hanning(N).^2) };
ih_L = fade_samples(ih_L,[50 50], fade_windows);
ih_R = fade_samples(ih_R,[50 50], fade_windows);

SH_ambisonic_binaural_decoder_DFE = zeros(length(SH_ambisonic_binaural_decoder(:,1,1))...
    ,length(SH_ambisonic_binaural_decoder(1,:,1)),length(SH_ambisonic_binaural_decoder(1,1,:)));
for i = 1:length(SH_ambisonic_binaural_decoder(:,1,1))
    SH_ambisonic_binaural_decoderEQ1(i,:,1) = conv(ih_L,SH_ambisonic_binaural_decoder(i,:,1));
    SH_ambisonic_binaural_decoderEQ1(i,:,2) = conv(ih_R,SH_ambisonic_binaural_decoder(i,:,2));
    
    SH_ambisonic_binaural_decoder_DFE(i,:,:) = SH_ambisonic_binaural_decoderEQ1(i...
        ,(length(ih_L)/2)+1:(length(SH_ambisonic_binaural_decoderEQ1(1,:,1)) - ...
        (length(ih_L)/2))+1,:); % THE +1 ensures a delay is not added of one sample (ttm 13.6.2018)
end

VL_hrirs_DFE = zeros(length(VL_hrirs(:,1,1))...
    ,length(VL_hrirs(1,:,1)),length(VL_hrirs(1,1,:)));
for i = 1:length(VL_hrirs(1,1,:))
    VL_hrirs_DFE1(1,:,i) = conv(ih_L,VL_hrirs(1,:,i));
    VL_hrirs_DFE1(2,:,i) = conv(ih_R,VL_hrirs(2,:,i));
    
    VL_hrirs_DFE(:,:,i) = VL_hrirs_DFE1(:...
        ,(length(ih_L)/2)+1:(length(VL_hrirs_DFE1(1,:,1)) - ...
        (length(ih_L)/2))+1,i); % THE +1 ensures a delay is not added of one sample (ttm 13.6.2018)
end

%% Plot results?
if dfe_plot_flag == 1
    comp_L = conv(ih_L, df_avg_L);
    
    plotColour = get(gca,'colororder'); figure;
    colour_frequency_plot(df_avg_L,Fs, window_size*2,  '-',plotColour(1,:), 1.8,[70 20000]);
    hold on
    colour_frequency_plot(ih_L,Fs, window_size*2, '-.',plotColour(2,:), 1.8,[70 20000]);
    colour_frequency_plot(comp_L,Fs, window_size*2, '-',[0 0 0], 1.4,[70 20000]);
    legend({'Diffuse-Field Response' ,...%(L)', 'Diffuse Field Response (R)', ...
        'Inverse Filter',...% (L)', 'Inverse Filter (R)',   ...
        'Result'}, 'location', 'southwest','FontSize',14);
    
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    set(gcf, 'Color', 'w');
    set(gca, 'fontsize', 14);
    xlim([70 20000]);
    ylim([-11 10]);
    disp('Ambisonic DFE complete!')
end
end

function [ azi,ele ] = tdesign2azi_ele( noOfPoints )
% tdesign2azi_ele()
%
% Converts spherical T-design quadrature to a matrix of azimuth and elevation.
% Azimuth angles 0 to 360 (going counterclockwise)
% Elevation angles -90 to 90 (going upwards)

file_ID = fopen(strcat('des',num2str(noOfPoints),'.txt'),'r'); %from http://neilsloane.com/sphdesigns/dim3/ (renamed text files)
raw_data = fscanf(file_ID,'%f');
xyz = reshape(raw_data,3,[])';

% initiate output arrays
azi = zeros(length(xyz),1); ele = zeros(length(xyz),1);
for i = 1:length(xyz)
    [azi(i),ele(i)] = xyz_to_tp_With_P_origin_at_0(xyz(i,1),xyz(i,2),xyz(i,3));
    if ele(i) == 90
        azi(i) = 0;
    elseif ele(i) == -90
        azi(i) = 0;
    end
end
end

