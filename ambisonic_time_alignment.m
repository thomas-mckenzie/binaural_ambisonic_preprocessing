function [VL_hrirs_TA] = ambisonic_time_alignment(VL_hrirs,Fs,ambisonic_order)
%{
Function to time-align high frequencies of virtual loudspeaker HRIRs for
Ambisonic reproduction.

Thomas McKenzie, University of York, 2019.

References:
See papers by Zaunschirm et al. (2018), as well as Evans et al. (1998),
Richter (2014).
%}

%%
% generate low-pass and high-pass crossover filters ----- IF lower than
% 4th order, only want to cross over to time aligned at ~2.5khz (so use crossover filter for 4th order).
if ambisonic_order < 4
    [filtLo,filtHi] = ambisonic_crossover(4,Fs);
else
    [filtLo,filtHi] = ambisonic_crossover(ambisonic_order,Fs);
end

% remove ITDs / time align
VL_hrirs_no_ITD = remove_ITD_HRIRs(VL_hrirs);

% crossover
for i = 1 : length(VL_hrirs(1,1,:))
    for j = 1:length(VL_hrirs(:,1,1))
        outputFilt_low(j,:,i) = filter(filtLo,1,VL_hrirs(j,:,i));
        outputFilt_High(j,:,i) = filter(filtHi,1,VL_hrirs_no_ITD(j,:,i));
    end
end

% Add the low freq and high freq together
VL_hrirs_TA = (outputFilt_low + outputFilt_High);

% filter delay compensation
gd = mean(grpdelay(filtLo));
for i = 1 : length(VL_hrirs(:,1,1))
    for j = 1:length(VL_hrirs(1,1,:))
        VL_hrirs_TA(i,:,j) = delayseq(VL_hrirs_TA(i,:,j)', -gd);
    end
end

disp('Ambisonic TA complete!');

end

function [VL_hrirs_ITD_removed] = remove_ITD_HRIRs(VL_hrirs)
VL_hrirs_ITD_removed = zeros(length(VL_hrirs(:,1,1)),length(VL_hrirs(1,:,1)),length(VL_hrirs(1,1,:)));

FcNorm = 500/(48000/2); % convert cutoff freq to normalised freq
XoverOrder = 8;
filtLo = fir1(XoverOrder,FcNorm,'low',hamming((XoverOrder+1)), 'noscale');

% remove ITD of each HRIR compared to the first
for i = 1 : length(VL_hrirs(1,1,:))
    for j = 1:length(VL_hrirs(:,1,1))
        hrir_1_f = filter(filtLo,1,VL_hrirs(1,:,1));
        hrir_i_f = filter(filtLo,1,VL_hrirs(j,:,i));
        hrirDelay = finddelay(hrir_1_f,hrir_i_f);

        if hrirDelay>=0
            VL_hrirs_ITD_removed(j,:,i) = [VL_hrirs(j,1+hrirDelay:length(VL_hrirs(j,:,i)),i), zeros(1,hrirDelay,1)];
        elseif hrirDelay < 0
            VL_hrirs_ITD_removed(j,:,i) = [zeros(1,abs(hrirDelay),1),VL_hrirs(j,1:length(VL_hrirs(j,:,i))-abs(hrirDelay),i)];
        end
    end
end

% Then remove all shared ITD:
inputHRIR1 = mean(mean(VL_hrirs,1),3);
inputHRIR2 = mean(mean(VL_hrirs_ITD_removed,1),3);
overallDelay = finddelay(inputHRIR1,inputHRIR2);

if overallDelay>=0
    VL_hrirs_ITD_removed = [VL_hrirs_ITD_removed(:,1+overallDelay:length(VL_hrirs_ITD_removed(1,:,1)),:), zeros(length(VL_hrirs(:,1,1)),overallDelay,length(VL_hrirs(1,1,:)))];
elseif overallDelay < 0
    VL_hrirs_ITD_removed = [zeros(length(VL_hrirs(:,1,1)),abs(overallDelay),length(VL_hrirs(1,1,:))),VL_hrirs_ITD_removed(:,1:length(VL_hrirs_ITD_removed(1,:,1))-abs(overallDelay),:)];
end
end


