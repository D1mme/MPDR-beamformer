%Author:    Dimme de Groot
%Date:      May 2024
%Descr:     This object implements the MPDR beamformer 

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

        %Some methods
        function window = methodSqrthann(obj)
            %See, i.e., Smith, J.O. Spectral Audio Signal Processing, http://ccrma.stanford.edu/~jos/sasp/, online book, 2011 edition.
            n = (0:1:obj.N_t-1).'; %note: obj.N_t is the length of the window
            m = n-obj.N_t/2;
            window = sqrt(ones(obj.N_t, 1).*(0.5+0.5*cos(2*pi/obj.N_t * m)));
        end
        function k_ax = methodK_ax(obj)
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
        function methodComputeBeamformerWeights(obj) %Eq. 3.23 of Brandstein 2001.       
            %Compute beamformed weights
            for k=1:obj.N_k
                R = obj.Rpsd(:,:,k) + obj.Rpsd_reg;
                lookVec = obj.look_vector(:,k);
                RinvV = R\lookVec;
                obj.w(:,k) = RinvV/(lookVec'*RinvV);
                obj.w(:,k) = obj.w(:,k);
            end
        end

        %Functions which are used during usage of beamformer
        function updateLookDirection(obj, look_location, receiver_location)
            %A function which computes the look direction vector obj.look_vector.
            %Inputs:    (1) look_location of size 1 x 3.
            %           (2) receiver_location of size Nr x 3, with Nr the number of receivers
            %Outputs:   an updated value of MPDRbeamformer.look_vector

            r = vecnorm(look_location-receiver_location,2,2);
            obj.look_vector = exp(-1j*obj.k_ax.*r)./(4*pi*r); %Note: no need to multiply by a constant or so since we are interested in a location
            %obj.look_vector = look_vector./vecnorm(look_vector); %Make the length equal to 1 per frequency bin: we now have a direction
        end

        function outputFrame = computeOutputFrame(obj, inputFrame, flag_bias_correction, VAD)
            %A function which computes the beamformed output frame. 
            %Inputs:    (1) the inputFrame - an unwindowned frame of signals as received at the mic of size window_length N_t x number of receivers N_r 
            %Outputs:   (1) output frame after beamforming
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
            obj.Rpsd = [];
            obj.Rpsd_reg = [];
        end
    end
end
