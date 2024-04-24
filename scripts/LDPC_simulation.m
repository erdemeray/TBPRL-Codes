%This script simulates the performance of TBP-RL LDPC (Type-Based Protograph Raptor-Like Low-Density Parity-Check) code over an AWGN (Additive White Gaussian Noise) channel.
% The LDPC code is defined by a Parity Check Matrix (PCM) loaded from a file.
% The simulation is performed for a range of SNR (Signal-to-Noise Ratio) values.
% The performance is evaluated in terms of Frame Error Rate (FER) and Bit Error Rate (BER).
clearvars

%% Set Parameters
NoI = 250; % number of decoding iterations
frames = 256; % number of erroneous frames
FER_sim_limit = 1e-2; %minimum FER limit
parallel_factor = 64;


rate_sim = [0.01,0.02,0.04]; % code rate - optimized for the rates between 0.2 and 0.01
%Calculate the SNR range for different rates SNR = (2^(Coderate*2)-1) 
SNR_start = floor(10*log10(2.^(rate_sim*2)-1)*10)/10; %in dB
SNR_NofE = 150;
stepsize = 0.05;
SNR_end = SNR_start + stepsize*SNR_NofE;
esno_dB = {SNR_start(1):stepsize:SNR_end(1);SNR_start(2):stepsize:SNR_end(2);SNR_start(3):stepsize:SNR_end(3)}; % Simulation SNR range 

    
    
% simulation result vectors
FER=zeros(numel(rate_sim),numel(esno_dB(1)));
BER=zeros(numel(rate_sim),numel(esno_dB(1)));

% Set the path of the PCM file.
file_name = 'Code/data/H_Rate_Adaptive_S1_50_S2_100_Expanded.mat';
addpath('util/','../data/','results/')

% initialize parallel cores
p=gcp('nocreate');
    
if isempty(p)
    num_of_cores = feature('numcores');
    parpool(num_of_cores);
else
    num_of_cores = p.NumWorkers;
end
%% Define the variables to be used

% check whether parallel_factor is integer multiple of erroneous frames
if ~(floor(frames/parallel_factor)==frames/parallel_factor)
    error("Number of erroneous frames must be an integer multiple of parallel factor.")
end

for j=1:numel(rate_sim)

    %channel parameters
    sigma = sqrt(10.^(-esno_dB{j}/10));
    Lc = 2*10.^(esno_dB{j}/10);

    H = get_H(rate_sim(j),file_name); 
    
    n = size(H,2);
    m=size(H,1);
    cfgLDPCDec = ldpcDecoderConfig(logical(H)); % decoder object
    
    %encode
    binary_CW = zeros(n,parallel_factor);
    x = -2 * binary_CW + 1; % transmitted codeword
    
    disp(['Used PCM : ' file_name]);
    
    disp(['Code Rate : ' num2str(rate_sim(j))]);
    
    disp(['Number of erroneous frames: ' num2str(frames)]);
    
    
    for i=1:numel(esno_dB{j})
        errors = 0;
        errors_frame = 0;
        total_frames=0;
        str_length = 0;
        while errors_frame<frames
            parfor frame = 1:(frames/parallel_factor)
                y = x + sigma(i)*randn(n,parallel_factor); % AWGN channel
    
                % calculate LLRs
                L = Lc(i)*y;
    
                %decode
                xh = ldpcDecode(L,cfgLDPCDec,NoI,'OutputFormat','whole','DecisionType','soft')<=0;
                
                %count the errors
                errors = errors + sum(xh ~= binary_CW,"all");
                if any(xh ~= binary_CW,'all')
                    errors_frame=errors_frame + sum(sum(xh ~= binary_CW,1)>0);
                end
            end
    
            total_frames=total_frames+frames;
            if total_frames > frames/ FER_sim_limit
                break;
            end
    
        end
    
        BER(j,i) = errors / total_frames / n;
        FER(j,i) = errors_frame / total_frames;
    
        str = sprintf('SNR: %.2f dB BER: %.5f FER: %.5f Erroneous / Total Frames: %i / %i',esno_dB{j}(i),BER(j,i),FER(j,i),errors_frame,total_frames);
        disp(str);
    
        if FER(j,i) < FER_sim_limit   
            disp(['Exiting the simulation due to the FER being less than the threshold ' num2str(FER_sim_limit) ' at the SNR ' num2str(esno_dB{j}(i))] );
            break;
        end
    end
end

%% save the results

save('results/results', '-regexp', '^(?!(H|x)$).')