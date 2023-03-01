function SH_ambisonic_binaural_decoder = generate_binaural_Ambisonic_decoder(ambisonic_order, VL_hrirs,loudspeaker_directions,dualband_flag)
%{
Function to generate a binaural Ambisonic decoder. Can be single band
(with basic SH channel weightings) or dual-band (with basic SH weightings
at low frequencies and Max rE weightings at high frequencies).

Thomas McKenzie, University of York, 2019.

References:
See Thomas McKenzie PhD thesis and Archontis Politis PhD thesis.

VL_hrirs should be in format:
[amount of loudspeakers x hrir x left&right]
%}

%%
if length(VL_hrirs(:,1,1)) < length(VL_hrirs(1,1,:))
    VL_hrirs = permute(VL_hrirs,[3 2 1]);
end

% get Fs
[~, Fs] = audioread(char(strcat('hrirs/azi_0_ele_0_DFC.wav')));

% % % just for the virtual loudspeakers:
% ambiDecoder from PolArch Ambisonic Library
ambiDeRV = ambiDecoder(loudspeaker_directions,'mmd',0,ambisonic_order);
ambiDeRE = ambiDecoder(loudspeaker_directions,'mmd',1,ambisonic_order);

shHrirRV = zeros(length(ambiDeRV(1,:)),length(VL_hrirs(1,:,1)),length(VL_hrirs(1,1,:)));
shHrirRE = zeros(length(ambiDeRE(1,:)),length(VL_hrirs(1,:,1)),length(VL_hrirs(1,1,:)));

for i = 1:length(VL_hrirs(1,1,:))
    for j = 1:length(VL_hrirs(:,1,1))
        shHrirRV(:,:,i) = shHrirRV(:,:,i) + (ambiDeRV(j,:)' * VL_hrirs(j,:,i));
        shHrirRE(:,:,i) = shHrirRE(:,:,i) + (ambiDeRE(j,:)' * VL_hrirs(j,:,i));
    end
end

%% dual band crossover
[filtLo,filtHi] = ambisonic_crossover(ambisonic_order,Fs);
gd = mean(grpdelay(filtLo));

% Filter the averages (RV --> low pass, Re --> high pass)
for i = 1:length(shHrirRV(:,1,1))
    for j = 1: length(shHrirRV(1,1,:))
        output_RVFilt(i,:,j) = filter(filtLo,1,shHrirRV(i,:,j));
        output_REFilt(i,:,j) = filter(filtHi,1,shHrirRE(i,:,j));
        shHrirDB(i,:,j) = (output_RVFilt(i,:,j) + output_REFilt(i,:,j));
        shHrirDB(i,:,j) = delayseq(shHrirDB(i,:,j)', -gd); % delay / advance to get rid of group delay from crossover
    end
end

if dualband_flag == 1
    SH_ambisonic_binaural_decoder = shHrirDB; % dual band
    %     disp('Dual-band binaural Ambisonic decoder generated!');
else
    SH_ambisonic_binaural_decoder = shHrirRV; % basic
    %     disp('Basic weighted binaural Ambisonic decoder generated!');
end
end

