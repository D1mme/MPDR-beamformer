%Author:    Dimme de Groot
%Date:      May 2024
%Descr:     Example 1 of how to use the MPDR beamformer object. In this example, there is only one interferer. 
%               This example validates the performance of the MPDR beamformer. The room impulse response (RIR) is the direct path (i.e. anechoic).
%               The look direction of the beamformer is the actual FFT(RIR), i.e. a perfect match.
%               In case you want to see how the MPDR beamformer performs with a near-perfect estimate of the noise PSD matrix, you can set VAD flag to true in line 71. 
%
%Note:      The code supports multiple interferers, just add more interferer locations. However, I didn't write code for multiple target locations. 

clear all 
close all

listen_flag = true;             %[-], set to true if you want to listen to the result
flag_bias_correction = true;    %[-], bias correction maskes the updates of Rpsd quicker for the first few frames. If you turn it off, the regulararisation dominates the first few frames, which is probably also fine

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User settings: MPDR settings and scenario: locations of target, receiver, interferer and received signal %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MPDR settings
fs = 16000;                     %[Hz], the sampling frequency
c = 342;                        %[m/s], the speed of sound
window_length =  0.016;         %[s], the window length of the analysis and synthesis window
pad_length = 0.016;             %[s], the amount of padding added
lambda = 0.98;                  %[-], used for updating the covariance  matrix: Rpsd_new = lambda*Rpsd_old + (1-lambda)*(xx^H), with x the observation vector    
regFactor = 1e-3;               %[W/hz or something], regularisation 
flag_full_axis = true;          %[-], either uses half of the frequency axis or the full frequency axis
analysis_window = "sqrthann";   %[-], the analysis window
synthesis_window = "sqrthann";  %[-], the synthesis window

%Scenario
locInterferer = [2.3,5.1,1];
locReceiver = [2.7, 4, 0.85;
               2.8, 4, 0.9;
               2.9, 4, 0.95;
               3.0, 4, 1.0;
               3.1, 4, 1.05;
               3.2, 4, 1.1;
               3.3, 4, 1.15;
               3.4, 4, 1.2];
locReceiver = locReceiver;
locTarget = [4, 5, 1];
plotLayout(locInterferer, locTarget, locReceiver)
[~, NN] = min(vecnorm(locReceiver-locTarget,2,2));
disp("The microphone which is nearest to the target is microphone number " + num2str(NN));
[audioRec, audioRecTAR, audioRecINT, audioPlayClean] = fnc_computeReceivedAudioAnechoic(locInterferer, locReceiver, locTarget,c , fs, window_length+pad_length);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise MPDR beamformer %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MPDR = MPDRbeamformer(c, fs, window_length, pad_length, lambda, regFactor, flag_full_axis, analysis_window, synthesis_window);

%Define look direction 
look_location = locTarget;          %The direction in which the beamformer "looks", this is defined with respect to....
receiver_location = locReceiver;    %... the receiver location 
MPDR.updateLookDirection(look_location, receiver_location);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform beamforming as if the frames come in on a frame-by-frame basis %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialise
stop_flag = false;
l = 0;
audioOut = zeros(size(audioRec,1),1);

%Loop over frames
while ~stop_flag

    %(1) Select a single frame; if an error is thrown, we are out of frames! Catch this and set stop flag to true
    try
        audioRecFrame = audioRec(l*MPDR.N_hop+1:l*MPDR.N_hop+MPDR.N_t,:);
    catch ME
        stop_flag = true;
    end

    %(2) If a frame is available the stop_flag is false.
    if ~stop_flag

        % For testing it is interesting to have a perfect voice activity detector (VAD). 
        % This essentially makes the MPDR an MVDR. To change from MPDR to MVDR set VAD=true in the if case below.
        % Note that this assumes a clean reference to be available!!!!
        if norm(audioRecTAR(l*MPDR.N_hop+1:l*MPDR.N_hop+MPDR.N_t, NN))>0.01
            VAD = false;  %To make MVDR: set this VAD to true
        else
            VAD = false; %Keep this one as is!
        end

        % Compute output frame and append to output signal
        audioOutFrame = MPDR.computeOutputFrame(audioRecFrame, flag_bias_correction, VAD); 
        audioOut(l*MPDR.N_hop+1:l*MPDR.N_hop+MPDR.N_t) = audioOut(l*MPDR.N_hop+1:l*MPDR.N_hop+MPDR.N_t) + audioOutFrame;
    
        l = l+1;
    end
end

%%%%%%%%%%%%%%%%
% Plot results %
%%%%%%%%%%%%%%%%
figure
plot(audioRec(:,NN))
hold on
plot(audioOut)
legend("Audio at microphone " + num2str((NN)), "Audio at output MPDR")
xlabel('Sample index')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute STOI at output of MPDR and at output of nearest neighbouring microphone %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
audioOut = audioOut(1:length(audioPlayClean));
[audioPlayCleanAligned, audioOutAligned] = alignsignals(audioPlayClean, audioOut, truncate=true);
disp("the estimated intelligibility (from 0 to 1, higher is better) is " + num2str(stoi(audioOutAligned, audioPlayCleanAligned, MPDR.fs)))

audioNN = audioRec(1:length(audioPlayClean), NN);
[audioPlayCleanAligned, audioNNAligned] = alignsignals(audioPlayClean, audioNN, truncate=true);
disp("the estimated intelligibility (from 0 to 1, higher is better) at the nearest neighbour is " + num2str(stoi(audioNNAligned, audioPlayCleanAligned, MPDR.fs)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Listen to audio signals %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if listen_flag
    soundsc(audioPlayCleanAligned, MPDR.fs);
    pause(length(audioPlayCleanAligned)/MPDR.fs)
    soundsc(audioNNAligned, MPDR.fs);
    pause(length(audioPlayCleanAligned)/MPDR.fs)
    soundsc(audioOutAligned, MPDR.fs);
end

%%%%%%%%%%%%%%%%%%
% Some functions %
%%%%%%%%%%%%%%%%%%
function plotLayout(locInterferer, locTarget, locReceiver)
    figure
    hold on
    scatter3(locInterferer(:,1), locInterferer(:,2), locInterferer(:,3), 'filled');
    scatter3(locTarget(:,1), locTarget(:,2), locTarget(:,3), 'filled');
    scatter3(locReceiver(:,1), locReceiver(:,2), locReceiver(:,3), 'filled');
    axis tight
    axis equal
    grid on
    legend('Interferer', 'Target', 'Microphones')
end

function [audioRec, audioRecTAR, audioRecINT, audioTar] = fnc_computeReceivedAudioAnechoic(locInterferer, locReceiver, locTarget, c, Fs, T)
    if nargin == 3
        c = 342;        %[m/s], sound velocity
        Fs = 16000;     %[Hz],  sampling frequency
        T = 0.256;      %[s],   estimated time-of-arrival plus bonus
    end

    Nrir = round(Fs*T);             %[-], number of samples in impulse response
    f_ax = 0:Fs/Nrir:Fs/2;          %[Hz], frequency axis
    f_ax = [f_ax(1:end-1),0, fliplr(-f_ax(2:end-1))];
    k_ax = 2*pi*f_ax/c;             %[?], wavenumber axis

    N_int = size(locInterferer,1);
    N_tar = size(locTarget, 1); %this should be one for the current implementation
    N_rec = size(locReceiver,1);

    audioTar = audioread('Data/sample_16kHz.wav');
    audioTar = [zeros(5*Fs,1); audioTar; zeros(3*Fs,1)];
    audioInt = 0.5*randn(length(audioTar), N_int);

    for i=1:N_rec
        for j=1:N_int
            r = norm(locInterferer(j,:)-locReceiver(i,:));
            RTF = 1/(4*pi*r)*exp(-1j.*k_ax*r);
            RIR = ifft(RTF);
            try
                audioRecINT(:,i) = audioRecINT(:,i) + conv(RIR,audioInt(:,j));
            catch 
                audioRecINT(:,i) = conv(RIR,audioInt(:,j));
            end
        end
        for j=1:N_tar
            r = norm(locTarget(j,:)-locReceiver(i,:));
            RTF = 1/(4*pi*r)*exp(-1j.*k_ax*r);
            RIR = ifft(RTF);
            try
                audioRecTAR(:,i) = audioRecTAR(:,i) + conv(RIR,audioTar(:,j));
            catch 
                audioRecTAR(:,i) = conv(RIR,audioTar(:,j));
            end
        end
    end
    audioRec = audioRecINT+audioRecTAR;

    %Give estimate SNR
    [~, NN] = min(vecnorm(locReceiver-locTarget,2,2));
    S = audioRecTAR(:,NN);
    N = audioRecINT(:,NN);
    INDX = find(abs(S)>0.01*max(abs(S)));
    S = S(INDX);
    N = N(INDX);
    SNR = 20*log10(norm(S)/norm(N));
    disp("estimated SNR at nearest neighbour is " + num2str(SNR) + " dB") 
end