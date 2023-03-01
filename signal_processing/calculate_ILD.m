function [ILD] = calculate_ILD(inputHRIR,Fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate ILD of an HRIR. 
% Works by high pass filtering at cutoff 2200hz (down to 1.5k), 
% then calculate magnitude response difference between left and right for 
% 30 different ERB bands, and getting the absolute mean of all differences 
% between 1.2k and 20k. 
% 
% ILD value in dB
% 
% Thomas McKenzie, University of York, 2018

if length(inputHRIR(:,1)) < length(inputHRIR(1,:)) 
    inputHRIR = inputHRIR';
end

hrirL = inputHRIR(:,1);
hrirR = inputHRIR(:,2);
if length(hrirL) < 1000
    Nfft = 4096;
else
Nfft = length(hrirL)*2;
end

lowFrequencyCutoff = 1200; %cutoff freq in hz
Fc = lowFrequencyCutoff/(Fs/2); % convert cutoff freq to normalised freq
XoverOrder = 128;
ripple = 50; % ripple in dB

filtHi = fir1(XoverOrder,Fc,'high',chebwin((XoverOrder+1),ripple), 'noscale');

gd = mean(grpdelay(filtHi));

hrirL_f = filter(filtHi,1,[hrirL; zeros(gd,1)]);
hrirR_f = filter(filtHi,1,[hrirR; zeros(gd,1)]);

hrirL_f = delayseq(hrirL_f, -gd); % delay / advance to get rid of group delay from crossover
hrirR_f = delayseq(hrirR_f, -gd); % delay / advance to get rid of group delay from crossover

% averaging over 30 ERB frequency bands (10 octaves so approx 1/3 octave bands). 
cfsAmount = 30;
cfs = (10.^((linspace((21.4*log10(4.37e-3*20+1)),...
    (21.4*log10(4.37e-3*20000+1)),cfsAmount))/21.4)-1)/4.37e-3; % from UNIVERSITY OF SURREY IOSR toolbox

HRTF_L = fft(hrirL_f, Nfft); % Get Fast Fourier transform
HRTF_R = fft(hrirR_f, Nfft); 

for i = 1:(cfsAmount-1)
frlowILD = round(cfs(i)*Nfft/Fs); % calculating ILD from LF cutoff to 20kHz
if frlowILD < 1
    frlowILD = 1;
end
frhigh = round(cfs(i+1)*Nfft/Fs);

magL1 = abs(HRTF_L(frlowILD:frhigh));
magR1 = abs(HRTF_R(frlowILD:frhigh));
ildPerBand(i) = mean(20*log10(magL1./magR1));

end

ILD = mean(ildPerBand(:));

end

