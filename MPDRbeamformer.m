%Author:    Dimme de Groot
%Date:      May 2024
%Descr:     This object implements the MPDR beamformer 
%
%           Note that the noise and or signal power spectral density matrices are iteratively estimated from the signals received at the microphones. 
%           Therefore, the MPDR beamformer is not necessarily expected to perform well for non-stationary noises (i.e. cars passing by, people talking, etc.)
%
%Methods:   
%       Create the beamformer object:           
%           MPDR = MPDRbeamformer(c, fs, window_length, pad_length, lambda, regFactor, flag_full_axis, analysis_window, synthesis_window);   
%               with:   c                   speed of sound [m/s], 
%                       fs                  sampling frequency [Hz]
%                       window_length       the length of the input frame in [s], 
%                       pad_length          the length of the padding in [s],
%                       lambda              a parameter used in updating the PSD estimates. Typically lambda in [0.97, 1)
%                       regFactor           the scaling of the regularisation identity matrix, 
%                       flag_full_axis      specifies whether or not to use the full frequency axis or only half
%                       analysis_window     a string specifying the analysis window (currently only: "sqrthann") 
%                       synthesis_window    a string specifying the synthesis window (currently only: "sqrthann") 
%           
%       MPDR beamformer: compute the output audio as obtained by the MPDR beamformer given an input audio signal
%           audioOut = MPDR_beamformer(obj, receiver_location, look_location, audioIn, resetPSD_flag, flag_bias_correction)
%               with:   look_location a 1 x 3 vector defining the target location in [m, m, m](i.e. where the desired source is located)
%                       receiver_location       a Nr x 3 vector defining the microphone locations in [m, m, m]
%                       audioIn                 the LENGTH INPUT x Nr matrix containing the audio received at the microphones
%                       resetPSD_flag           an optional input argument (default: false) to reset the PSD matrices (i.e. reset the MPDR object)
%                       flag_bias_correction    an optional input argument (default: true) to perform bias correction when updating the PSD matrix estimates. 
%               note:   in principle, we are only interested in the distances from the microphones to the look_location. 
%                       Thus, there is some redundant information in specifying the locations; I did so anyway as for my use case this was the easiest. 
%                       It is easy to modify the function to take distances directly: see the method updateLookDirection(look_location, receiver_location)        
%                      
%       MVDR beamformer: compute the output audio as obtained by an MVDR beamformer with oracle voice activity detection given an input audio signal and the clean input audio signal. 
%           audioOut = MVDR_beamformer(obj, receiver_location, look_location, audioIn, audioInClean, resetPSD_flag, flag_bias_correction)
%               with:   look_location           a 1 x 3 vector defining the target location in [m, m, m](i.e. where the desired source is located)
%                       receiver_location       a Nr x 3 vector defining the microphone locations in [m, m, m]
%                       audioIn                 the LENGTH INPUT x Nr matrix containing the audio received at the microphones
%                       audioInClean            the LENGTH INPUT x Nr matrix containing the clean(!) audio received at the microphones. Clean means that there is no noise: only the audio produced by the desired speaker.
%                       resetPSD_flag           an optional input argument (default: false) to reset the PSD matrices (i.e. reset the MPDR object)
%                       flag_bias_correction    an optional input argument (default: true) to perform bias correction when updating the PSD matrix estimates. 
%               Note:   In practical cases we of course dont have the clean input signal available. This function mostly shows how to use the VAD parameter   

classdef MPDRbeamformer < handle
    properties
        c                   %[m/s], speed of sound
        fs                  %[Hz],  sampling frequency

        window_length_act   %[s], actual window length
        N_t                 %[-], window length in samples      
        pad_length_act      %[s], actual padding length
        N_pad               %[-], padding length in samples
        N_hop               %[-], hop length in samples        
        w_analysis          %[-], analysis window
        w_synthesis         %[-], synthesis window  

        lambda              %[-], Parameter of exponential window which controls the update of the covariance matrix.
        regFactor           %[?], Parameter of regularisation applied to estimated covariance matrics
       
        Rpsd                %[?], The estimated covariance matrix (excluding regularisation)
        Rpsd_reg            %[?], The regularisation term for the estimated covariance matrix    

        look_vector         %[-], the look vector. Essentially specifies which direction we like to keep.

        flag_full_axis      %[-], True for full frequency axis [-Fs/2, Fs/2). False for [0, Fs/2]
        k_ax                %[?], wavenumber
        N_k                 %[-], length of k_axis

        %For debugging
        w                   %[-], the beamformer weights
        FFTinputFrame       %[-]      
        FFToutputFrame      %[-]

        Rpsd_bias_corr      %In case the bias correction flag is true, updating Rpsd also requires tracking the "biased" version of the Rpsd. In case the flag is false, Rpsd = Rpsd_bias_corr
        frame_counter       %[-], counts the number of times Rpsd has been updated
    end
    methods
        function obj = MPDRbeamformer(c, fs, window_length, pad_length, lambda, regFactor, flag_full_axis, analysis_window, synthesis_window)
            if nargin == 0
                c = 342;
                fs = 16000;

                window_length = 0.016;
                pad_length = 0.016;
                    
                lambda = 0.9;
                regFactor = 1e-5;

                flag_full_axis = false;      
                analysis_window = "sqrthann";
                synthesis_window = "sqrthann";
            end
            obj.c = c;
            obj.fs = fs;
            
            obj.N_t = 2^nextpow2(floor(fs.*window_length));    
            if pad_length ~= 0
                obj.N_pad = 2^nextpow2(floor(fs.*pad_length));   
            else
                obj.N_pad = 0;
            end
            obj.window_length_act = obj.N_t/fs;
            obj.pad_length_act = obj.N_pad/fs;
            
            obj.lambda = lambda;
            obj.regFactor = regFactor;
            obj.flag_full_axis = flag_full_axis;    

            if analysis_window == "sqrthann"
               obj.N_hop = obj.N_t/2;
               obj.w_analysis = methodSqrthann(obj); 
            else
                disp("currently only sqrthann is supported as window :S")
            end

            if synthesis_window == "sqrthann"
               obj.w_synthesis = methodSqrthann(obj); 
            else
                disp("currently only sqrthann is supported as window :S")
            end

            obj.k_ax = methodK_ax(obj);
            obj.N_k = length(obj.k_ax);
        end

        %Some methods used internally
        function window = methodSqrthann(obj)
            % This function computes thes sqaure root hanning window
            % See, i.e., Smith, J.O. Spectral Audio Signal Processing, http://ccrma.stanford.edu/~jos/sasp/, online book, 2011 edition.
            n = (0:1:obj.N_t-1).'; %note: obj.N_t is the length of the window
            m = n-obj.N_t/2;
            window = sqrt(ones(obj.N_t, 1).*(0.5+0.5*cos(2*pi/obj.N_t * m)));
        end
        function k_ax = methodK_ax(obj)
            % this function computes the wavenumber axis
            N_k = obj.N_t + obj.N_pad;                              %[-], total number of frequency bins (even)  
            if obj.flag_full_axis
                k_ax = (-N_k/2:1:N_k/2-1)/N_k*2*pi*obj.fs/obj.c; 
                k_ax = fftshift(k_ax); 
            else
                k_ax = (0:1:N_k/2)/N_k*2*pi*obj.fs/obj.c;          %[rad/m], wave number -> only half of the axis [0, Fs/2] is relevant
            end
        end
        function methodUpdateRpsd(obj, FFTinputFrame, flag_bias_correction)
            
            if nargin == 2
                flag_bias_correction = false;
            end
            
            if flag_bias_correction
                if isempty(obj.frame_counter)
                    obj.frame_counter = 1;
                else
                    obj.frame_counter = obj.frame_counter + 1;
                end
            end

            %Loop over k-axis
            for k=1:obj.N_k
                audioRecFrame_k = FFTinputFrame(k,:).'; 
                Rpsd_new = (1-obj.lambda)*(audioRecFrame_k*audioRecFrame_k') + obj.lambda*obj.Rpsd_bias_corr(:,:,k);
                obj.Rpsd_bias_corr(:,:,k) = Rpsd_new;   
                
                if flag_bias_correction
                    obj.Rpsd(:,:,k) = Rpsd_new/(1-obj.lambda^obj.frame_counter);
                else
                    obj.Rpsd(:,:,k) = Rpsd_new; 
                end  
            end
        end
        function methodComputeBeamformerWeights(obj)
            % this function computes the beamformer weights, see Eq. (3.23) of Brandstein & Ward, microphone arrays, 2001.  
            for k=1:obj.N_k
                R = obj.Rpsd(:,:,k) + obj.Rpsd_reg;
                lookVec = obj.look_vector(:,k);
                RinvV = R\lookVec;
                obj.w(:,k) = RinvV/(lookVec'*RinvV);
                obj.w(:,k) = obj.w(:,k);
            end
        end

        function updateLookDirection(obj, look_location, receiver_location)
            %A function which computes the look direction vector obj.look_vector.
            %Inputs:    (1) look_location of size 1 x 3.
            %           (2) receiver_location of size Nr x 3, with Nr the number of receivers
            %Outputs:   an updated value of MPDRbeamformer.look_vector

            r = vecnorm(look_location-receiver_location,2,2);
            obj.look_vector = exp(-1j*obj.k_ax.*r)./(4*pi*r); %Note: no need to multiply by a constant or so since we are interested in a location
        end

        function outputFrame = computeOutputFrame(obj, inputFrame, flag_bias_correction, VAD)
            %A function which computes the beamformed output frame. 
            % Inputs:    (1) the inputFrame - an unwindowned frame of signals as received at the mic of size window_length N_t x number of receivers N_r 
            % Outputs:   (1) output frame after beamforming
            if nargin == 2
                flag_bias_correction = false;
                VAD = false;
            end
            if nargin == 3 
                VAD = false;
            end

            Nr = size(inputFrame,2);                            %number of receivers
            inputFrame = obj.w_analysis.*inputFrame;            %window the input frame

            %Initialise PSD matrices if needed
            if isempty(obj.Rpsd_reg)
                obj.Rpsd_reg = eye(Nr)*obj.regFactor;
            end
            if isempty(obj.Rpsd)
                obj.Rpsd = zeros(Nr, Nr, obj.N_k); 
                obj.Rpsd_bias_corr = zeros(Nr, Nr, obj.N_k);
            end

            %Take fft and select correct part of axis
            obj.FFTinputFrame = fft(inputFrame, obj.N_t+obj.N_pad, 1); %Take columnwise fft and zero pad to proper length
            if ~obj.flag_full_axis
                obj.FFTinputFrame = obj.FFTinputFrame(1:obj.N_k,:);
            end

            %Update PSD matrix
            if ~VAD
                obj.methodUpdateRpsd(obj.FFTinputFrame, flag_bias_correction);
            end

            %Compute beamformer weights
            obj.methodComputeBeamformerWeights;

            %Compute output frame
            obj.FFToutputFrame = obj.w'.*obj.FFTinputFrame;
            if obj.flag_full_axis           
                outputFrame = ifft(sum(obj.FFToutputFrame,2), 'symmetric');
            else
                outputFrame = sum(ifft([obj.FFToutputFrame; zeros(obj.N_t+obj.N_pad-obj.N_k,Nr)], 'symmetric'),2);
            end
            outputFrame = outputFrame(1:obj.N_t).*obj.w_synthesis;
        end

        function resetPSD(obj) 
            % This function resets the power spectral density matrix and is handy when using the object multiple times independently
            obj.Rpsd = [];
            obj.Rpsd_reg = [];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Main function: MPDR beamformer %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function audioOut = MPDR_beamformer(obj, receiver_location, look_location, audioIn, resetPSD_flag, flag_bias_correction)
            %This function takes the receiver (microphone) location, the targhet source location, and the input audio.
            %Optionally, it also takes resetPSD_flag (default: false) for resetting the PSD matrix estimates and flag_bias_correction (default: true)
                            
            if nargin == 4
                resetPSD_flag = false;
                flag_bias_correction = true;
            end
            if nargin == 5
                flag_bias_correction = true;
            end
        
            %Reset PSD matrices
            if resetPSD_flag    
                obj.resetPSD;
            end

            %Define look direction  
            obj.updateLookDirection(look_location, receiver_location);
        
            %Initialise
            stop_flag = false;                      % Becomes true when we are out of frames
            l = 0;                                  % keeps track of the frame number        
            audioOut = zeros(size(audioIn,1),1);    % Placeholder for the output audio
            
            %Note: for the MPDR, VAD is always false. 
            %   In case you want to use an MVDR beamformer, VAD should be true for input frames which have voice activity.
            %   For synthetic data, this is easily accomplished by detecting if the clean (noiseless) signal has energy for a given frame.  
            VAD = false;    
        
            %Loop over frames
            while ~stop_flag
                try %If we are out of input frames, an error is caught and the stop_flag is set to true.
                    audioRecFrame = audioIn(l*obj.N_hop+1:l*obj.N_hop+obj.N_t,:);
                catch ME
                    stop_flag = true;
                end
            
                if ~stop_flag
                    %Compute output frame and append to outut signal
                    audioOutFrame = obj.computeOutputFrame(audioRecFrame, flag_bias_correction, VAD); 
                    audioOut(l*obj.N_hop+1:l*obj.N_hop+obj.N_t) = audioOut(l*obj.N_hop+1:l*obj.N_hop+obj.N_t) + audioOutFrame;
                
                    l = l+1;
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Main function: MVDR beamformer %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function audioOut = MVDR_beamformer(obj, receiver_location, look_location, audioIn, audioCleanIn, resetPSD_flag, flag_bias_correction)
            %This function takes the receiver (microphone) location, the targhet source location, and the input audio.
            %Optionally, it also takes resetPSD_flag (default: false) for resetting the PSD matrix estimates and flag_bias_correction (default: true)
                            
            if nargin == 5
                resetPSD_flag = false;
                flag_bias_correction = true;
            end
            if nargin == 6
                flag_bias_correction = true;
            end
        
            %Reset PSD matrices
            if resetPSD_flag    
                obj.resetPSD;
            end

            %Define look direction  
            obj.updateLookDirection(look_location, receiver_location);
        
            %Initialise
            stop_flag = false;                      % Becomes true when we are out of frames
            l = 0;                                  % keeps track of the frame number        
            audioOut = zeros(size(audioIn,1),1);    % Placeholder for the output audio
            
            %Note: for the MPDR, VAD is always false. 
            %   In case you want to use an MVDR beamformer, VAD should be true for input frames which have voice activity.
            %   For synthetic data, this is easily accomplished by detecting if the clean (noiseless) signal has energy for a given frame.  
            VAD = false;    
        
            %Loop over frames
            while ~stop_flag
                try %If we are out of input frames, an error is caught and the stop_flag is set to true.
                    audioRecFrame = audioIn(l*obj.N_hop+1:l*obj.N_hop+obj.N_t,:);
                    audioRecCleanFrame = audioCleanIn(l*obj.N_hop+1:l*obj.N_hop+obj.N_t,:);
                catch ME
                    stop_flag = true;
                end
            
                if ~stop_flag
                    %Oracle voice activity detection
                    if norm(audioRecCleanFrame) ~= 0
                        VAD = true;
                    else
                        VAD = false;
                    end

                    %Compute output frame and append to outut signal
                    audioOutFrame = obj.computeOutputFrame(audioRecFrame, flag_bias_correction, VAD); 
                    audioOut(l*obj.N_hop+1:l*obj.N_hop+obj.N_t) = audioOut(l*obj.N_hop+1:l*obj.N_hop+obj.N_t) + audioOutFrame;
                
                    l = l+1;
                end
            end
        end
    end
end
