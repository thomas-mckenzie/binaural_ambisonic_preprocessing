function db_values = colour_frequency_plot(input_signal,Fs, window_size, plot_style, plot_colour, line_width,x_limits)
% Thomas McKenzie, University of York, 2019. 

input_signal_FD = fft(input_signal, window_size); % obtain Frequency Domain version of input signal
f_low = round(x_limits(1)*window_size/Fs); % Compute freq bins for x-axis limits
f_high = round(x_limits(2)*window_size/Fs);
freq_vector = Fs/window_size:Fs/window_size:Fs; % Frequency vector for plotting

db_values = 20*log10(abs(input_signal_FD(f_low:f_high)));
semilogx(freq_vector(f_low:f_high),db_values, plot_style,'Color',plot_colour, 'linewidth', line_width); % plot
grid on;
