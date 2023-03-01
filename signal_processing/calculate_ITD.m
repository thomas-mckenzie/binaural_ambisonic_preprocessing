function [ITD,ITD_samples] = calculate_ITD(inputHRIR,Fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate ITD of an HRIR
% Works by using finddelay of a 1500 Hz low-passed version of 
% the HRIR. 
% 
% ITD value in miliseconds (ms)
% 
% Thomas McKenzie, University of York, 2018

if length(inputHRIR(:,1)) < length(inputHRIR(1,:)) 
    inputHRIR = inputHRIR';
end

hrirL = inputHRIR(:,1);
hrirR = inputHRIR(:,2);

Fc = 1500/(Fs/2); % convert cutoff freq to normalised freq
XoverOrder = 32*16;

filtLo = fir1(XoverOrder,Fc,'low',hamming((XoverOrder+1)), 'noscale');

gd = mean(grpdelay(filtLo));

hrirL_f = filter(filtLo,1,[hrirL; zeros(gd,1)]);
hrirR_f = filter(filtLo,1,[hrirR; zeros(gd,1)]);

hrirL_f = delayseq(hrirL_f, -gd); % delay / advance to get rid of group delay from crossover
hrirR_f = delayseq(hrirR_f, -gd); % delay / advance to get rid of group delay from crossover


ITD_samples =  finddelay(hrirL_f,hrirR_f);
ITD = ITD_samples / Fs * 1000; 
end

