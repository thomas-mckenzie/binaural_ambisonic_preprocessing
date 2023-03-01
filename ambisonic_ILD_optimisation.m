function [SH_ambisonic_binaural_decoder_AIO,VL_hrirs_augmented_final,ILD_augmentation_gains,ILD_augmentation_gains_final] = ambisonic_ILD_optimisation(ambisonic_order,SH_ambisonic_binaural_decoder,input_VL_hrirs,reference_VL_hrirs,Fs,loudspeaker_directions, dualband_flag)
%{
Function to augment high frequencies of virtual loudspeaker HRIRs for
Ambisonic reproduction, in order to improve the Ambisonic reproduction
accuracy of interaural level differences.

Thomas McKenzie, University of York, 2019.

References:
"Interaural level difference optimisation of binaural Ambisonic rendering,"
Applied Sciences, vol. 9, no. 6, 2019, T. McKenzie, D. T. Murphy, and G.
Kearney.
%}

%%
maxIterations = 101;

% Initial renders
[ILD_augmentation_gains] = calculate_ILD_differences(reference_VL_hrirs,Fs,SH_ambisonic_binaural_decoder,loudspeaker_directions,ambisonic_order);
i = 1;

% Second render
ILD_augmentation_gains_product = prod(ILD_augmentation_gains,1);
VL_hrirs_augmented = augment_ILD(input_VL_hrirs,Fs,ambisonic_order,ILD_augmentation_gains_product);
shHrirDB_augmented = generate_binaural_Ambisonic_decoder(ambisonic_order, VL_hrirs_augmented,loudspeaker_directions,dualband_flag); %encode virtual loudspeakers to SH domain
[ILD_augmentation_gains(i+1,:)] = calculate_ILD_differences(reference_VL_hrirs,Fs,shHrirDB_augmented,loudspeaker_directions,ambisonic_order);
i = 2;

% While there is still change between the ILD stretch values, keep iterating
while (round(mean(ILD_augmentation_gains(i,:)),5) ~= round(mean(ILD_augmentation_gains(i-1,:)),5))
    ILD_augmentation_gains_product = prod(ILD_augmentation_gains,1);
    VL_hrirs_augmented = augment_ILD(input_VL_hrirs,Fs,ambisonic_order,ILD_augmentation_gains_product);
    shHrirDB_augmented = generate_binaural_Ambisonic_decoder(ambisonic_order, VL_hrirs_augmented,loudspeaker_directions,dualband_flag); %encode virtual loudspeakers to SH domain
    [ILD_augmentation_gains(i+1,:)] = calculate_ILD_differences(reference_VL_hrirs,Fs,shHrirDB_augmented,loudspeaker_directions,ambisonic_order);
    
    i = i+1;
    if i > maxIterations %if it's iterated 100 times then quit the while loop
        disp(strcat('iterations exceed maximum of ',num2str(maxIterations)));
        break
    end
end

% Once final values for augmentation have been found, make the new HRIRs:
ILD_augmentation_gains_final = prod(ILD_augmentation_gains);
VL_hrirs_augmented_final = augment_ILD(input_VL_hrirs,Fs,ambisonic_order,ILD_augmentation_gains_final);
SH_ambisonic_binaural_decoder_AIO = generate_binaural_Ambisonic_decoder(ambisonic_order, VL_hrirs_augmented_final,loudspeaker_directions,dualband_flag); %encode virtual loudspeakers to SH domain

disp(strcat('Ambisonic ILD Optimisation complete! Number of iterations = ',num2str(i)));
end

function ILD_augmentation_gains = calculate_ILD_differences(VL_hrirs,Fs,SH_ambisonic_binaural_decoder,loudspeaker_directions,ambisonic_order)

% Generate ambisonic HRIR for the position of each loudspeaker.
% Calculate ILD of the HRIR and the ambisonic HRIR
for i = 1:length(loudspeaker_directions)
    encodedAmbisonicGains = encodeHOA_N3D(ambisonic_order,1,[loudspeaker_directions(i,1) loudspeaker_directions(i,2)]);
    ambisonicHRIR(1,:) = encodedAmbisonicGains * SH_ambisonic_binaural_decoder(:,:,1);
    ambisonicHRIR(2,:) = encodedAmbisonicGains * SH_ambisonic_binaural_decoder(:,:,2);
    
    ild_HRIR(i) = calculate_ILD(VL_hrirs(:,:,i),Fs);
    ild_Ambi(i) = calculate_ILD(ambisonicHRIR,Fs);
end

% calculate difference in absolute ILD (in dB)
ild_difference = abs(ild_HRIR) - abs(ild_Ambi);

% calculate gain factor corresponding to difference
ILD_augmentation_gains = (10.^(ild_difference/20));

for i = 1 : length(VL_hrirs(1,1,:))
    if loudspeaker_directions(i,1) == 0 || loudspeaker_directions(i,1) == 180 % Don't change median plane HRIRs
        ILD_augmentation_gains(i) = 1;
    end
end
end

function [VL_hrirs_ILD_augmented_dualband] = augment_ILD(VL_hrirs,Fs,ambisonic_order,gainValues)

for i = 1:length(VL_hrirs(1,1,:))
    ildH1(i) = calculate_ILD(VL_hrirs(:,:,i),Fs);
end

for i = 1 : length(VL_hrirs(1,1,:))
    if gainValues(i) == 1 % if HRIR should stay the same
        VL_hrirs_augmented(1,:,i) = VL_hrirs(1,:,i);
        VL_hrirs_augmented(2,:,i) = VL_hrirs(2,:,i);
    else
        if ildH1(i) > 0 % if ILD needs to get augmented on left side
            VL_hrirs_augmented(1,:,i) = VL_hrirs(1,:,i);
            VL_hrirs_augmented(2,:,i) = VL_hrirs(2,:,i)/ gainValues(i);
        elseif ildH1(i) < 0 % if ILD needs to get augmented on right side
            VL_hrirs_augmented(2,:,i) = VL_hrirs(2,:,i);
            VL_hrirs_augmented(1,:,i) = VL_hrirs(1,:,i)/ gainValues(i);
        end
    end
end

%% Normalise augmented ILD HRIRs to same RMS as original HRIRs

% Correcting gains so that stretched IRs are same RMS level as unstretched IRs
% Note: commenting this out can create much bigger ILDs across the circle,
% BUT this then makes PSD worse
for i = 1 : length(VL_hrirs(:,1,1))
    for j = 1 : length(VL_hrirs(1,1,:))
        rmsH1(i,j) = rms(VL_hrirs(i,:,j)); %%% NORMALISE TO RMS of whole HRIR
        rmsS1(i,j) = rms(VL_hrirs_augmented(i,:,j));
    end
end

meanRmsH = mean(rmsH1);
meanRmsS = mean(rmsS1);

for i = 1:length(VL_hrirs(1,1,:))
    VL_hrirs_augmented(:,:,i) =  VL_hrirs_augmented(:,:,i) * meanRmsH(i)/ meanRmsS(i)   ;
end

% generate low-pass and high-pass crossover filters
lfCutoffOrder = 3;
if ambisonic_order < lfCutoffOrder
    [filtLo,filtHi] = ambisonic_crossover(1500,Fs); % if ambi order is lower than 3, then only implement aio at frequencies above 1500hz. This value is taken due to Middlebrooks1991
else
    [filtLo,filtHi] = ambisonic_crossover(ambisonic_order,Fs);
end

% filter
for i = 1 : length(VL_hrirs(1,1,:))
    for j = 1:length(VL_hrirs(:,1,1))
        outputFilt_low(j,:,i) = filter(filtLo,1,VL_hrirs(j,:,i));
        outputFilt_High(j,:,i) = filter(filtHi,1,VL_hrirs_augmented(j,:,i));
    end
end

% Add the low freq and high freq together
VL_hrirs_ILD_augmented_dualband = (outputFilt_low + outputFilt_High);

% compensate for delay caused by filters
gd = mean(grpdelay(filtLo));
for i = 1 : length(VL_hrirs(:,1,1))
    for j = 1:length(VL_hrirs(1,1,:))
        VL_hrirs_ILD_augmented_dualband(i,:,j) = delayseq(VL_hrirs_ILD_augmented_dualband(i,:,j)', -gd);
    end
end
end
