function [loudspeaker_directions, VL_hrirs,Fs] = load_hrirs_lebedev(ambisonic_order)
%{
Function to load virtual loudspeaker HRIRs for Ambisonic reproduction.
These are Bernschutz wavs (diffuse-field equalised by Tom McKenzie 2018)

Thomas McKenzie, University of York, 2019.

References:
See papers by Bernschutz (2013) for HRIR measurement methods,
Thomas McKenzie PhD thesis for HRIR diffuse-field equalisation method.
%}

%%
switch ambisonic_order
    case 1
        disp('First Order: Forward facing Octahedron');
        loudspeaker_directions = [0 0; 180 0; 270 45; 270 -45; 90 45; 90 -45]; %forward facing octa
    case 2
        disp('Second Order: 14 pt. Lebedev Grid');
        loudspeaker_directions =      [0,180,90,270,0,0,45,135,315,225,45,135,315,225]'; %14 pt. Lebedev Grid
        loudspeaker_directions(:,2) = [0,0,0,0,-90,90,-35,-35,-35,-35,35,35,35,35]';
    case 3
        disp ('Third Order: 26 pt. Lebedev Grid');
        loudspeaker_directions = [0 0 90 180 270 45	135	225	315 0 45 90	135	180 225	270	315	45 135	225	315	0 90 180 270 0]';
        loudspeaker_directions(:,2) = [90 45 45 45 45 35 35 35 35 0 0 0 0 0 0 0 0 -35 -35 -35 -35 -45 -45 -45 -45 -90]';
    case 4
        disp ('Fourth Order: 38 pt. Lebedev Grid'); % note - the \hrirs folder is missing some necessary measurements (as the 38 pt. grid is not nested inside the 50 pt. grid) so revert to 50 pt. grid instead.
%         loudspeaker_directions =      [0,180,90,270,0,0,45,135,315,225,45,135,315,225,62,118,298,242,28,152,332,208,0,180,0,180,0,180,0,180,90,270,90,270,90,270,90,270]';
%         loudspeaker_directions(:,2) = [0,0,0,0,-90,90,-35,-35,-35,-35,35,35,35,35,0,0,0,0,0,0,0,0,-62,-62,62,62,-28,-28,28,28,-62,-62,62,62,-28,-28,28,28]';
        disp ('Fourth Order: 38 pt. Lebedev Grid unavailable so using 50 pt. Lebedev Grid instead');
        loudspeaker_directions =       [0  45	135	225 315	0  90   180 270 45 135	225 315	18  72 108  162 198 252 288 342 0	45  90 135 180 225  270 315 18  72  108	162	198 252	288	342 45  135 225 315 0   90	180	270	 45 135 225 315  0]';
        loudspeaker_directions(:,2) = [90	64	64	64	64  45 45	45	45	35 35	35	35	18	18 18	18	18	18	18  18	0	0	0	0	0	0	0    0	-18	-18	-18 -18	-18	-18	-18 -18	-35	-35	-35 -35	-45	-45	-45 -45	-64	-64	-64 -64	-90]';
    case 5
        disp ('Fifth Order: 50 pt. Lebedev Grid');
        loudspeaker_directions =       [0  45	135	225 315	0  90   180 270 45 135	225 315	18  72 108  162 198 252 288 342 0	45  90 135 180 225  270 315 18  72  108	162	198 252	288	342 45  135 225 315 0   90	180	270	 45 135 225 315  0]';
        loudspeaker_directions(:,2) = [90	64	64	64	64  45 45	45	45	35 35	35	35	18	18 18	18	18	18	18  18	0	0	0	0	0	0	0    0	-18	-18	-18 -18	-18	-18	-18 -18	-35	-35	-35 -35	-45	-45	-45 -45	-64	-64	-64 -64	-90]';
    otherwise
        error('Unsupported Ambisonic order');
end

[testHRIR,Fs] = audioread('hrirs/azi_0_ele_0_DFC.wav');
winSize = length(testHRIR(:,1));
VL_hrirs = zeros(2,winSize+1,length(loudspeaker_directions));

for i = 1:length(loudspeaker_directions)
    % If it can't find exact number, first tries the closest multiple of 2 for azi then
    % closest multiple of 2 for ele then just goes closest multiple of 2 for
    % both if can't find either
    if exist(strcat('hrirs/azi_',num2str(loudspeaker_directions(i,1)),'_ele_',num2str(loudspeaker_directions(i,2)),'_DFC.wav')) == 2
        hrirs = (audioread(char(strcat('hrirs/azi_',num2str(loudspeaker_directions(i,1)),'_ele_',num2str(loudspeaker_directions(i,2)),'_DFC.wav'))))';
    elseif exist(strcat('hrirs/azi_',num2str((round(loudspeaker_directions(i,1)/2)*2)),'_ele_',num2str(loudspeaker_directions(i,2)),'_DFC.wav')) == 2
        loudspeaker_directions(i,1) = (round(loudspeaker_directions(i,1)/2)*2);
        hrirs = (audioread(char(strcat('hrirs/azi_',num2str(loudspeaker_directions(i,1)),'_ele_',num2str(loudspeaker_directions(i,2)),'_DFC.wav'))))';
    elseif exist(strcat('hrirs/azi_',num2str(loudspeaker_directions(i,1)),'_ele_',num2str(round(loudspeaker_directions(i,2)/2)*2),'_DFC.wav')) == 2
        loudspeaker_directions(i,2) = (round(loudspeaker_directions(i,2)/2)*2);
        hrirs = (audioread(char(strcat('hrirs/azi_',num2str(loudspeaker_directions(i,1)),'_ele_',num2str(loudspeaker_directions(i,2)),'_DFC.wav'))))';
    elseif exist(strcat('hrirs/azi_',num2str(loudspeaker_directions(i,1)),'_ele_',num2str(loudspeaker_directions(i,2)),'_DFC.wav')) == 0
        dhalf = loudspeaker_directions(i,:);
        dhalf = (round(dhalf/2)*2);
        if dhalf(1,1) == 360
            dhalf(1,1) = 0;
        end
        loudspeaker_directions(i,:) = dhalf;
        hrirs = (audioread(char(strcat('hrirs/azi_',num2str(loudspeaker_directions(i,1)),'_ele_',num2str(loudspeaker_directions(i,2)),'_DFC.wav'))))';
    end
    
    % Place in matrix VL_hrirs
    VL_hrirs(:,1:winSize,i) = hrirs;
end
end