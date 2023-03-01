%{
This script allows the testing of binaural Ambisonic decoders. It plots
the diffuse-field response of the decoders, and then compares the
binaural Ambisonic rendering versus the standard (non-Ambisonic) HRIRs in
folder \hrirs.
This is in three ways: perceptual spectral difference (PSD), interaural
level difference (ILD) and interaural time difference (ITD).

Run this script after running load_ambisonic_configuration.m

Thomas McKenzie, University of York, 2019.

References:
See Thomas McKenzie PhD thesis.
%}

%% Load in HRIRs

winSize = length(SH_ambisonic_binaural_decoder(1,:,1));

% Create a list of recorded files in directory
s = dir(strcat('hrirs/*.wav')); % s is structure array
file_list = {s.name}'; % convert the name field into cell array of strings.
num_meas = length(s); % Number of measurements

ambiMatrix_NPP = zeros(winSize,num_meas,2); % initialise matrices
ambiMatrix = zeros(winSize,num_meas,2); % initialise matrices
hrirMatrix = zeros(winSize,num_meas, 2); % Create empty array for HRIRs
azimuth_H = zeros(1, num_meas);
elevation_H = zeros(1, num_meas);

for i = 1:num_meas
    % Put each pair of HRIRs into a big array
    [hrirMatrix(1:winSize-1,i,:)] = audioread(char(strcat('hrirs/', file_list(i,:))));
    
    % Extract Azimuth and Elevation angle based on filename
    filenamestr = char(file_list(i,:)); % Get current filename
    IndexAzi = strfind(file_list(i,:), 'azi_'); % Find the text 'azi_'
    azimuth_H(i) = sscanf(filenamestr(1,cell2mat(IndexAzi) + ...
        length('azi_'):end), '%g', 1); % Get azimuth value
    IndexEle = strfind(file_list(i,:), 'ele_'); % Find the text 'ele_'
    elevation_H(i) = sscanf(filenamestr(1,cell2mat(IndexEle) + ...
        length('ele_'):end), '%g', 1); % Get elevation value
    
    % generate ambisonic HRIRs
    encodedAmbisonicGains = encodeHOA_N3D(ambisonic_order,1,[azimuth_H(i) elevation_H(i)]);
    
    ambisonicHRIR_NPP(1,:) = encodedAmbisonicGains * SH_ambisonic_binaural_decoder_NPP(:,:,1);
    ambisonicHRIR_NPP(2,:) = encodedAmbisonicGains * SH_ambisonic_binaural_decoder_NPP(:,:,2);
    
    ambisonicHRIR(1,:) = encodedAmbisonicGains * SH_ambisonic_binaural_decoder(:,:,1);
    ambisonicHRIR(2,:) = encodedAmbisonicGains * SH_ambisonic_binaural_decoder(:,:,2);
    
    ambiMatrix_NPP(:,i,:) = ambisonicHRIR_NPP';
    ambiMatrix(:,i,:) = ambisonicHRIR';
end

Nfft = length(hrirMatrix(:,1,1));

%% view diffuse-field responses
[~,~,df_avg_L_NPP,~] = ambisonic_diffuse_field_equalisation(ambisonic_order,SH_ambisonic_binaural_decoder_NPP,VL_hrirs_NPP,Fs,0);
[~,~,df_avg_L,~] = ambisonic_diffuse_field_equalisation(ambisonic_order,SH_ambisonic_binaural_decoder,VL_hrirs,Fs,0);

plotColour = get(gca,'colororder');

figure;
set(gcf, 'Position',  [100, 10, 760, 800]); 
hold on
subplot(4,2,1);
dB_val_NPP = colour_frequency_plot(df_avg_L_NPP,Fs,Nfft*2,'-',plotColour(1,:),1.8,[70 20000]);
xlabel('Frequency (Hz)');ylabel('Magnitude (dB)');xlim([70 20000]);
title('Diffuse-field response (NPP)');
ylim([-20 10]);

subplot(4,2,2);
dB_val = colour_frequency_plot(df_avg_L,Fs, Nfft*2,  '-',plotColour(2,:), 1.8,[70 20000]);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); xlim([70 20000]);
title('Diffuse-field response');
ylim([-20 10]);
ylim([min([dB_val_NPP;dB_val])-2 max([dB_val_NPP;dB_val])+2]);

subplot(4,2,1);
ylim([min([dB_val_NPP;dB_val])-2 max([dB_val_NPP;dB_val])+2]);

%% Perceptual Spectral Difference

solid_angle_weighting = GetVoronoiPlotandSolidAng(azimuth_H, elevation_H, 0);
solid_angle_weighting = solid_angle_weighting/(sum(solid_angle_weighting));

f.fs = 48000; f.nfft = length(hrirMatrix(:,1,1)); f.minFreq = 20; f.maxFreq = 20000;

% using the mckenzie2022 predicted binaural colouration model -- see
% https://www.amtoolbox.org/models.php and 
% T. McKenzie, C. Armstrong, L. Ward, D. Murphy, and G. Kearney.
% Predicting the colouration between binaural signals. Appl. Sci.,
% 12(2441), 2022.
[~,PavgSpecDiff_NPP] = mckenzie2022(ambiMatrix_NPP,hrirMatrix,0,f,[],solid_angle_weighting,0,0.01);
[~,PavgSpecDiff] = mckenzie2022(ambiMatrix,hrirMatrix,0,f,[],solid_angle_weighting,0,0.01);

for i = 1: num_meas
    specDiffAmbiSAW_NPP(i,:) = PavgSpecDiff_NPP(:,i,:) * solid_angle_weighting(i);
    specDiffAmbiSAW(i,:) = PavgSpecDiff(:,i,:) * solid_angle_weighting(i);
end

weightedAmbiSpecDiffL_NPP = sum(specDiffAmbiSAW_NPP(:,1));
weightedAmbiSpecDiffR_NPP = sum(specDiffAmbiSAW_NPP(:,2));
weightedAmbiSpecDiffLandR_NPP = ((weightedAmbiSpecDiffL_NPP+weightedAmbiSpecDiffR_NPP)/2);
disp(strcat('Average PSD NPP= ',num2str(weightedAmbiSpecDiffLandR_NPP)))

weightedAmbiSpecDiffL = sum(specDiffAmbiSAW(:,1));
weightedAmbiSpecDiffR = sum(specDiffAmbiSAW(:,2));
weightedAmbiSpecDiffLandR = ((weightedAmbiSpecDiffL+weightedAmbiSpecDiffR)/2);
disp(strcat('Average PSD = ',num2str(weightedAmbiSpecDiffLandR)))

specDiffAmbi_NPP = squeeze(PavgSpecDiff_NPP);
specDiffAmbi = squeeze(PavgSpecDiff);

x = azimuth_H; y = elevation_H;
for i = 1: length(x)
    if x(i) > 180
        x(i) = x(i)-360;
    end
end

z = mean(specDiffAmbi_NPP,2)';
z2 = mean(specDiffAmbi,2)';

subplot(4,2,3);
heatmap_plot(x,y,z);%caxis([1 2.5]);
title('PSD (NPP)');
c2 = colorbar;
c2.Label.String = 'PSD (sones)';
ylabel ('Elevation (°)')
caxis([min([z z2]) max([z z2])]);

subplot(4,2,4);
heatmap_plot(x,y,z2);%caxis([1 2.5]);
title('PSD');
c2 = colorbar;
c2.Label.String = 'PSD (sones)';
ylabel ('Elevation (°)')
caxis([min([z z2]) max([z z2])]);

%% ILD

for i = 1:num_meas
    ILD_hrir(i) = calculate_ILD(squeeze(hrirMatrix(:,i,:)),Fs);
    ILD_ambi(i) = calculate_ILD(squeeze(ambiMatrix(:,i,:)),Fs);
    ILD_ambi_NPP(i) = calculate_ILD(squeeze(ambiMatrix_NPP(:,i,:)),Fs);
end

ILD_D_ambi_NPP = abs(ILD_hrir - ILD_ambi_NPP);% to plot just absolute difference
ILD_D_ambi = abs(ILD_hrir - ILD_ambi);% to plot just absolute difference

mILD_D_ambi_NPP = (abs(ILD_hrir - ILD_ambi_NPP));
mILD_D_ambi = (abs(ILD_hrir - ILD_ambi));
for i = 1: num_meas
    ILD_D_ambi_SAW(1,i) = mILD_D_ambi(1,i) * solid_angle_weighting(i);
    ILD_D_ambi_SAW_NPP(1,i) = mILD_D_ambi_NPP(1,i) * solid_angle_weighting(i);
end

meanILD_D_ambiSAW_NPP = sum(ILD_D_ambi_SAW_NPP);
disp(strcat('Average \DeltaILD NPP = ',num2str(meanILD_D_ambiSAW_NPP)))
meanILD_D_ambiSAW = sum(ILD_D_ambi_SAW);
disp(strcat('Average \DeltaILD = ',num2str(meanILD_D_ambiSAW)))

subplot(4,2,5);
heatmap_plot(x,y,ILD_D_ambi_NPP);colorbar; % caxis([0 7]);
title('ILD (NPP)');
c2 = colorbar;
c2.Label.String = '\DeltaILD (dB)';
caxis([min([ILD_D_ambi_NPP ILD_D_ambi]) max([ILD_D_ambi_NPP ILD_D_ambi])]);

subplot(4,2,6);
heatmap_plot(x,y,ILD_D_ambi);colorbar; % caxis([0 7]);
title('ILD');
c2 = colorbar;
c2.Label.String = '\DeltaILD (dB)';
caxis([min([ILD_D_ambi_NPP ILD_D_ambi]) max([ILD_D_ambi_NPP ILD_D_ambi])]);

%% ITD

for i = 1:num_meas
    ITD_hrir(i) = calculate_ITD(squeeze(hrirMatrix(:,i,:)),Fs);
    ITD_ambi(i) = calculate_ITD(squeeze(ambiMatrix(:,i,:)),Fs);
    ITD_ambi_NPP(i) = calculate_ITD(squeeze(ambiMatrix_NPP(:,i,:)),Fs);
end

ITD_D_ambi_NPP = abs(ITD_hrir - ITD_ambi_NPP); % to see just pos values
ITD_D_ambi = abs(ITD_hrir - ITD_ambi); % to see just pos values

mITD_D_ambi_NPP = (abs(ITD_hrir - ITD_ambi_NPP));
mITD_D_ambi = (abs(ITD_hrir - ITD_ambi));
for i = 1: num_meas
    ITD_D_ambi_SAW_NPP(1,i) = mITD_D_ambi_NPP(1,i) * solid_angle_weighting(i);
    ITD_D_ambi_SAW(1,i) = mITD_D_ambi(1,i) * solid_angle_weighting(i);
end

meanITD_D_ambiSAW_NPP = sum(ITD_D_ambi_SAW_NPP);
disp(strcat('Average \DeltaITD NPP= ',num2str(meanITD_D_ambiSAW_NPP)))
meanITD_D_ambiSAW = sum(ITD_D_ambi_SAW);
disp(strcat('Average \DeltaITD = ',num2str(meanITD_D_ambiSAW)))

subplot(4,2,7);
heatmap_plot(x,y,ITD_D_ambi_NPP);colorbar; % caxis([0 0.4]);
title('ITD (NPP)');
c2 = colorbar;
c2.Label.String = '\DeltaITD (ms)';
caxis([min([ITD_D_ambi_NPP ITD_D_ambi]) max([ITD_D_ambi_NPP ITD_D_ambi])]);

subplot(4,2,8);
heatmap_plot(x,y,ITD_D_ambi);colorbar; % caxis([0 0.4]);
title('ITD');
c2 = colorbar;
c2.Label.String = '\DeltaITD (ms)';
caxis([min([ITD_D_ambi_NPP ITD_D_ambi]) max([ITD_D_ambi_NPP ITD_D_ambi])]);
