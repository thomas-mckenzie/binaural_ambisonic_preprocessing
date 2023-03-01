%{
This test script illustrates how to generate binaural Ambisonic renders
using an excerpt of the EigenScape freely available fourth-order
Ambisonic database: https://zenodo.org/record/1012809
Green, Marc, and Damian Murphy. "EigenScape: A database of spatial
acoustic scene recordings." Applied Sciences 7.11 (2017): 1204

Run this script after running load_ambisonic_configuration.m

Thomas McKenzie, University of York, 2019.

References:
See Thomas McKenzie PhD thesis.
%}

%%
% ambisonic_soundscape = audioread('eigenscape_beach.wav');
ambisonic_soundscape = audioread('eigenscape_train.wav');

% convert from SN3D normalisation to N3D
ambisonic_soundscape = convert_N3D_SN3D(ambisonic_soundscape, 'sn2n');

% initialise matrices
binaural_ambisonic_render_NPP = zeros(length(SH_ambisonic_binaural_decoder(1,:,1))+length(ambisonic_soundscape)-1,length(SH_ambisonic_binaural_decoder(1,1,:)));
binaural_ambisonic_render = zeros(length(SH_ambisonic_binaural_decoder(1,:,1))+length(ambisonic_soundscape)-1,length(SH_ambisonic_binaural_decoder(1,1,:)));

% convolve each channel of the encoded signal with the decoder signal and sum the result
for i = 1:length(SH_ambisonic_binaural_decoder(:,1,1))
    binaural_ambisonic_render_NPP(:,1) = binaural_ambisonic_render_NPP(:,1) +  conv(SH_ambisonic_binaural_decoder_NPP(i,:,1),ambisonic_soundscape(:,i) );
    binaural_ambisonic_render_NPP(:,2) = binaural_ambisonic_render_NPP(:,2) +  conv(SH_ambisonic_binaural_decoder_NPP(i,:,2),ambisonic_soundscape(:,i) );
    
    binaural_ambisonic_render(:,1) = binaural_ambisonic_render(:,1) +  conv(SH_ambisonic_binaural_decoder(i,:,1),ambisonic_soundscape(:,i) );
    binaural_ambisonic_render(:,2) = binaural_ambisonic_render(:,2) +  conv(SH_ambisonic_binaural_decoder(i,:,2),ambisonic_soundscape(:,i) );
end

% compare the two binaural decoders by listening to both binaural renders consecutively
soundsc([binaural_ambisonic_render_NPP; binaural_ambisonic_render],Fs);