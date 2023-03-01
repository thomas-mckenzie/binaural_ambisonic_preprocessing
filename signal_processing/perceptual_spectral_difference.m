function [diff, avgDiff_ERB, offset_F] = perceptual_spectral_difference(A, B, solidAngle, freq, limit, offset_F)
%{
Function to compare the spectra of two datasets of binaural signals using a
perceptually motivated calculation.

Calum Armstrong and Thomas McKenzie, University of York, 2018. 

References:
"A perceptual spectral difference model for binaural signals," in AES 145th
Convention, New York, USA, 2018, C. Armstrong, T. McKenzie, D. T. Murphy, 
and G. Kearney.
%}


%{
- Function operation - 

This compares the spectra of A and B in terms of PERCEPTUALLY WEIGHTED
error. For multi-dimensional inputs the Difference is made along the first
dimension. Averages are output for each column of data.

Perceptual error takes into account the lesser importance of quieter
sounds. Frequency bins are weighted with respect to the ISO 226 loudness
curves for an average listeneing level of 75 dB SPL. A comtribution half as
loud is deemed half as important.

The perceptual average difference is further weighted with respoect to ERB
bandwidth. Simply, this reduces the contibution to the average calculation
of higher frequency components where our ears are less sensitive. It's like
a logorythmic type average.

An iterative optimisation process is used to find the input normalisation
which results in the lowest error metric. This is generally somewhere
around the point at which the two input signals have the same mean value -
but can easily vary by a few dB. The optimum normalisation is different for
perceptual / absolute error metrics.

If you do not wish any normalisation to be applied, optional input offset_F
should be set to 'none'. A value of 0 will skip the optimal normalisation
process and just equate the mean values. Any other value will specify the
additional normalisation (in terms of dB gain) applied to B AFTER
normalisation with respect to the mean values of the two inputs.

*** A and B are expected to be input on a dB scale
%}


%% SETUP

% Check and format input vectors
if size(A, 1) == 1
    A = A';
end
if size(B, 1) == 1
    B = B';
end
if size(freq, 1) == 1
    freq = freq';
end
if size(solidAngle, 1) == 1
    solidAngle = solidAngle';
end

% Test for input errors
if size(A) ~= size(B)
    error('Input matricies A and B must be the same size');
end

% unless offset_F = 'none' normalise means
if not(exist('offset_F', 'var') && strcmp(offset_F, 'none'))
    % normalise mean values of A and B
    avgA = mean(A(:)); avgB = mean(B(:));
    absNormOff = avgA - avgB;
    B = B + absNormOff;
else
    offset_F = 0;
end

% Define what limits number of iterations
if limit >= 1
    maxIterations = limit;
    minResolution = 0;
    lastIteration = limit;
elseif limit > 0
    maxIterations = 1000;
    minResolution = limit;
else
    error('Positive numeric value expected for input variable "limit"');
end

% Set initial values for iterative normalization
initInc = 0.2;
scale = 0.4;
offset = zeros(maxIterations, 1);
avgPError = zeros(maxIterations, 1);

%% FUNCTION SETUP

% Useful variable for later
sizeA = size(A);

% ISO 226 Declarations
iso226SPL = zeros(30, 91);
iso226Freq = zeros(30, 91);
Y = zeros(length(freq), 91);
EL = zeros(length(freq), 91);

% For all ISO standardised listening levels (ISO 226: 0-90dB SPL)
for l = 1:91
    % Save equal loudness contour
    [iso226SPL(1:29, l), iso226Freq(1:29, l)] = iso226(l-1);
    iso226Freq(30, l) = 20000;
    iso226SPL(30, l) = iso226SPL(1, l);
    
    % Fit curve to equal loudness contour
    iso226Fit = fit(iso226Freq(:, l), iso226SPL(:, l), 'pchip');
    
    % Interpolate to input frequency bins and remove equivilant 1KHz
    % loudness offset
    Y(:, l) = iso226Fit(freq) - (l - 1);
    
    % Save the offset required in dB to equate the loudness of any
    % frequency bin to that of 1KHz for a given absolute loudness
    % (0-90dB SPL)
    % ... A.K.A. flip the equal loudness contour 1KHz offset!
    EL(:, l) = -Y(:, l) ;
end

% Calculate ERB bandwidths
ERB = 0.108.*freq + 24.7;

% Calculate ERB weights and repeat for input matrix multiplication
ERBWeights_temp = (1./ERB) ./ max(1./ERB);
ERBWeights = repmat(ERBWeights_temp, [1, sizeA(2:end)]);

%% ITERATE TO FIND PERCEPTUALLY OPTIMUM INPUT NORMALISATION

if nargin < 6
    
    %         update = figure(); hold on;
    
    [~, temp, jump] = get_spectral_difference(A, B, 0);
    temp_weighted = temp .* solidAngle';
    avgPError(1) = sum(temp_weighted(:))/2 ./ sum(solidAngle(:));
    
    %         Update figure
    %         figure(update);
    %         plot(0, avgPError(1), 'x','linewidth', 1.5,'MarkerSize',10); hold on;
    %         drawnow
    
    % For each iteration
    for i = 2:maxIterations
        switch i
            case 2 % 1st (No offset)
                
                offset(i) = jump;
                
            case 3 % 2nd (Try out positive offset)
                
                offset(i) = jump + initInc;
                
            otherwise % All others
                
                % If there was a reduction in perceptual error with
                % respect to the previous offset, carry on increasing
                % the offset in that direction (up or down)
                if avgPError(i-1) < avgPError(i-2)
                    
                    offset(i) = offset(i-1) + initInc;
                    
                    % If there was no improvement...
                else
                    % If this is only the third iteration, switch the
                    % offset direction, but keep the offset resolution
                    % the same
                    if i == 4
                        
                        initInc = initInc * -1;
                        offset(i) = jump + initInc;
                        
                        % Otherwise, carry on for two more samples then
                        % switch the offset direction and reduce the
                        % resultion (narrowing in on optimum value)
                    else
                        
                        initInc = initInc * -scale;
                        
                        % If the resolution has got too fine, quit
                        if abs(initInc) < minResolution
                            lastIteration = i-1;
                            break;
                        else
                            offset(i) = offset(i-1) + initInc;
                        end
                    end
                end
        end
        
        % Perform Difference with calculated offset and save result
        [~, temp, ~] = get_spectral_difference(A, B, offset(i));
        temp_weighted = temp .* solidAngle';
        avgPError(i) = sum(temp_weighted(:))/2 ./ sum(solidAngle(:));
        
        %             Update figure
        %             figure(update);
        %             plot(offset(i), avgPError(i), 'x','linewidth', 1.5,'MarkerSize',10); hold on;
        %             drawnow
        
    end
    
    % Extract unique results
    [offset_U, ia, ~] = unique(offset(1:lastIteration));
    avgPError_U = avgPError(ia);
    
    % Find approximate minimum
    [~, minIndex] = min(avgPError_U);
    
    offset_F = offset_U(minIndex);
    %         fprintf('Optimum Normalisation Found: %d\n', offset_F);
end


%% COMPARE SPECTRA WITH PERCEPTUALLY OPTIMAL NORMALISATION

% Perform Difference
[diff, avgDiff_ERB, ~] = get_spectral_difference(A, B, offset_F);

%     close(update);

%% GET SPECTRAL DIFFERENCE

    function [pDiff, avgPDiff_ERB, jump] = get_spectral_difference(A, B, offset_F)
        %% Re-Normalise inputs and amplify to 'Listening Level' amplitudes
        %  ... Average value set to 75dB SPL
        
        B = B + offset_F;
        
        % normalise to 75dB
        combined_matrices = ([A B]);
        meanValue = mean(combined_matrices(:));
        %                 meanValue = mean2(combined_matrices);
        A = A + (75-meanValue);
        B = B + (75-meanValue);
        
        %% Account for Equal Loudness
        
        % Save matricies that select the corect loudness contour offsets to
        % use for each frequency bin for each input spectrum. This is the
        % integer rounded values of the input matricies with a maximum
        % value of 90 and a minimum value of 0.
        LC_A = min(round(A), 90);
        LC_A = max(LC_A, 0);
        LC_B = min(round(B), 90);
        LC_B = max(LC_B, 0);
        
        EL_A = EL(sub2ind(size(EL), repmat((1:sizeA(1))', [1, sizeA(2:end)]), LC_A+1));
        EL_B = EL(sub2ind(size(EL), repmat((1:sizeA(1))', [1, sizeA(2:end)]), LC_B+1));
        
        % Account for equal loudness / convert to phones scale
        A_EL = A + EL_A;
        B_EL = B + EL_B;
        
        
        %% Convert to Sones
        A_EL_Sones = 2.^((A_EL-40)/10);
        B_EL_Sones = 2.^((B_EL-40)/10);
        
        %% FIND LOUDNESS WEIGHTED DIFFERENCES
        
        % Perceptual Difference
        pDiff = (B_EL_Sones - A_EL_Sones);
        avgPDiff_ERB = sum(ERBWeights .* abs(pDiff)) ./ sum(ERBWeights);
        erbSumA = sum(ERBWeights .* (A_EL_Sones)) ./ sum(ERBWeights);
        erbSumB = sum(ERBWeights .* (B_EL_Sones)) ./ sum(ERBWeights);
        %         jump = mean2(sum(ERBWeights .* (A_EL_Sones)) ./ sum(ERBWeights)) - mean2(sum(ERBWeights .* (B_EL_Sones)) ./ sum(ERBWeights));
        jump = mean(erbSumA(:)) - mean(erbSumB(:));
        
        
    end
end



function [spl, freq] = iso226(phon)
%
% Generates an Equal Loudness Contour as described in ISO 226
%
% Usage:  [SPL FREQ] = ISO226(PHON);
% 
%         PHON is the phon value in dB SPL that you want the equal
%           loudness curve to represent. (1phon = 1dB @ 1kHz)
%         SPL is the Sound Pressure Level amplitude returned for
%           each of the 29 frequencies evaluated by ISO226.
%         FREQ is the returned vector of frequencies that ISO226
%           evaluates to generate the contour.
%
% Desc:   This function will return the equal loudness contour for
%         your desired phon level.  The frequencies evaulated in this
%         function only span from 20Hz - 12.5kHz, and only 29 selective
%         frequencies are covered.  This is the limitation of the ISO
%         standard.
%
%         In addition the valid phon range should be 0 - 90 dB SPL.
%         Values outside this range do not have experimental values
%         and their contours should be treated as inaccurate.
%
%         If more samples are required you should be able to easily
%         interpolate these values using spline().
%
% Author: Jeff Tackett 03/01/05



%                /---------------------------------------\
%%%%%%%%%%%%%%%%%          TABLES FROM ISO226             %%%%%%%%%%%%%%%%%
%                \---------------------------------------/
f = [20 25 31.5 40 50 63 80 100 125 160 200 250 315 400 500 630 800 ...
     1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 10000 12500];

af = [0.532 0.506 0.480 0.455 0.432 0.409 0.387 0.367 0.349 0.330 0.315 ...
      0.301 0.288 0.276 0.267 0.259 0.253 0.250 0.246 0.244 0.243 0.243 ...
      0.243 0.242 0.242 0.245 0.254 0.271 0.301];

Lu = [-31.6 -27.2 -23.0 -19.1 -15.9 -13.0 -10.3 -8.1 -6.2 -4.5 -3.1 ...
       -2.0  -1.1  -0.4   0.0   0.3   0.5   0.0 -2.7 -4.1 -1.0  1.7 ...
        2.5   1.2  -2.1  -7.1 -11.2 -10.7  -3.1];

Tf = [ 78.5  68.7  59.5  51.1  44.0  37.5  31.5  26.5  22.1  17.9  14.4 ...
       11.4   8.6   6.2   4.4   3.0   2.2   2.4   3.5   1.7  -1.3  -4.2 ...
       -6.0  -5.4  -1.5   6.0  12.6  13.9  12.3];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%Error Trapping
if((phon < 0) | (phon > 90))
    disp('Phon value out of bounds!')
    spl = 0;
    freq = 0;
else
    %Setup user-defined values for equation
    Ln = phon;

    %Deriving sound pressure level from loudness level (iso226 sect 4.1)
    Af=4.47E-3 * (10.^(0.025*Ln) - 1.15) + (0.4*10.^(((Tf+Lu)/10)-9 )).^af;
    Lp=((10./af).*log10(Af)) - Lu + 94;

    %Return user data
    spl = Lp;  
    freq = f;
end

end