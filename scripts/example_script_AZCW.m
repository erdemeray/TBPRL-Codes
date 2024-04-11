% This script simulates the performance of TBP-RL LDPC (Type-Based Protograph Raptor-Like Low-Density Parity-Check) code over an AWGN (Additive White Gaussian Noise) channel.
% The LDPC code is defined by a Parity Check Matrix (PCM) loaded from a file.
% The simulation is performed for a range of SNR (Signal-to-Noise Ratio) values.
% The performance is evaluated in terms of Frame Error Rate (FER) and Bit Error Rate (BER).

%% Parameters to be set
NoI=500; % number of decoding iterations
frames=128; % number of erroneous frames
FER_sim_limit = 1e-2; % minimum FER limit

esno_dB = -7.7:0.05:-7.6; % Simulation SNR range
rate_sim = 0.2; % code rate - optimized for the rates between 0.2 and 0.01

%esno_dB = -11.25:0.05:-11; % Simulation SNR range
%rate_sim = 0.1; % code rate - optimized for the rates between 0.2 and 0.01

parallel_factor = 32;%number of frames decoded in parallel.

% Set the path of the PCM file.
file_name = '../PCM/H_AZCW.mat';

% Multithreading/GPU flag
% Note: In one decoder, "parallel_factor" number of frames are decoded in parallel. 
% If you want you can decode in multiple decoders for speed-up. 
% This is especially beneficial when you have multiple GPUs installed in the machine. 
multithread_flag = true; % use parallel decoders 
GPU_flag = false; % use the GPU decoder 

%% Define the variables to be used

p=gcp('nocreate');

num_threads = 1;

if multithread_flag ==true
    if GPU_flag ==true
        if isempty(p)
            num_threads = gpuDeviceCount("available");
            parpool(num_threads);
        else
            num_threads = p.NumWorkers;
        end
    else
        if isempty(p)
            num_threads = feature('numcores');
            parpool(num_threads);
        else
            num_threads = p.NumWorkers;
        end
    end
end

% check whether parallel_factor is integer multiple of erroneous frames
if ~(floor(frames/parallel_factor)==frames/parallel_factor)
    error("Number of erroneous frames must be an integer multiple of parallel factor.")
end

H = get_H(rate_sim,file_name);

n = size(H,2);
m=size(H,1);
cfgLDPCDec = ldpcDecoderConfig(logical(H)); % decoder object

%encode
binary_CW = zeros(n,parallel_factor);
x = -2 * binary_CW + 1; % transmitted codeword

%channel parameters
sigma = sqrt(0.5 * 10.^(-esno_dB/10)); % noise std.dev.
Lc = 4*10.^(esno_dB/10);

% simulation result vectors
FER=zeros(NoI,numel(esno_dB));
BER=zeros(NoI,numel(esno_dB));


for i=1:numel(esno_dB)
    errors = 0;
    errors_frame = 0;
    total_frames=0;
    str_length = 0;
    while errors_frame<frames
        for frame = 1:(frames/parallel_factor)

            if multithread_flag==true
                % Multithread - GPU or CPU
                frame_errors_thread=zeros(1,num_threads);
                bit_errors_thread=zeros(1,num_threads);

                parfor thread_idx=1:num_threads
                    y = x + sigma(i)*randn(n,parallel_factor); % AWGN channel

                    % calculate LLRs
                    L = Lc(i)*y;

                    if GPU_flag==true
                        L=gpuArray(L);
                    end

                    %decode
                    xh = ldpcDecode(L,cfgLDPCDec,NoI,'OutputFormat','whole','DecisionType','soft')<=0;

                    %count the errors
                    bit_errors_thread(thread_idx) = sum(xh ~= binary_CW,"all");
                    if any(xh ~= binary_CW,'all')
                        frame_errors_thread(thread_idx) = sum(sum(xh ~= binary_CW,1)>0);
                    end
                end

                errors = errors + sum(bit_errors_thread);
                errors_frame = errors_frame + sum(frame_errors_thread);


            else % Single Thread
                y = x + sigma(i)*randn(n,parallel_factor); % AWGN channel

                % calculate LLRs
                L = Lc(i)*y;

                if GPU_flag==true
                    L=gpuArray(L);
                end

                %decode
                xh = ldpcDecode(L,cfgLDPCDec,NoI,'OutputFormat','whole','DecisionType','soft')<=0;

                %count the errors
                errors = errors + sum(xh ~= binary_CW,"all");
                if any(xh ~= binary_CW,'all')
                    errors_frame = errors_frame + sum(sum(xh ~= binary_CW,1)>0);
                end
            end
        end

        total_frames=total_frames+frames*num_threads;

    end

    BER(i) = errors / total_frames / n;
    FER(i) = errors_frame / total_frames;

    str = sprintf('SNR: %.2f dB BER: %.5f FER: %.5f Erroneous / Total Frames: %i / %i',esno_dB(i),BER(i),FER(i),errors_frame,total_frames);
    disp(str);

    if ((FER(i) < FER_sim_limit)||(total_frames > frames/FER_sim_limit))
        disp(['Exiting the simulation due to the FER being less than the threshold ' num2str(FER_sim_limit) ' at the SNR ' num2str(esno_dB(i))] );
        break;
    end
end

delete(gcp('nocreate'));