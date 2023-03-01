function [matrix_output_fft, freq_vector_fft] = fft_matrix(matrix_input, Fs, Nfft, freq_range)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate the single sided frequency spectrum of two matrices
% of HRIRs for a specified frequency range. Returns FFT of input matrix as
% the absolute FFT in dB for the specified frequency range with the
% associated frequency vector. 
%
% ====================== Input Parameters =======================
%
% matrix_input:             input matrix of HRIRs - in format [HRIR, no of
%                           measurements, left and right]
%
% Fs:                       Sample rate
%
% Nfft:                     Length of FFT
%
% freq_range:               axis limits in Hz - in format [fr_low fr_high]
%
% Tom McKenzie - returns fft in dB of two matrices 2018

% Take FFT of matrices
fft_matrix_input = fft(matrix_input, Nfft); % Get Fast Fourier transform

% Compute freq bins for x-axis limits
fr_low = round(freq_range(1)*Nfft/Fs);
fr_high = round(freq_range(2)*Nfft/Fs);

% Get absolute values for frequency bins
fft_abs_matrix_input = abs(fft_matrix_input(fr_low:fr_high,:,:));

% Get values in dB
matrix_output_fft = 20*log10(fft_abs_matrix_input);

% Frequency vector for plotting
f = 0:Fs/Nfft:Fs-(Fs/Nfft);
freq_vector_fft = f(fr_low:fr_high);

end