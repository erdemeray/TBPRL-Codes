% This function loads a PCM (Parity Check Matrix) from a file and adjusts it based on the desired rate.
% The PCM is then converted to a sparse matrix for efficiency.
%
% @param rate_sim: The desired rate for the simulation.
% @param file_name: The path of the file containing the PCM.
% @return H: The adjusted and sparse PCM to be used for decoding.

function [H] = get_H(rate_sim, file_name)

    % Load the PCM from the specified file.
    load(file_name);

    % Calculate the lifting factor, which is a scaling factor used in protograph-based LDPC codes.
    lifting_factor = 5000;

    % Calculate the number of variable nodes (VNs) in the PCM based on the desired rate and the lifting factor.
    num_of_PCM_VNs = floor((4*lifting_factor)/rate_sim);

    % Calculate the number of check nodes (CNs) in the PCM based on the number of VNs and the lifting factor.
    num_of_PCM_CNs = num_of_PCM_VNs - 4 *  lifting_factor;

    % Adjust the PCM to have the calculated number of CNs and VNs.
    H = H(1:num_of_PCM_CNs,1:num_of_PCM_VNs);

    % Convert the PCM to a sparse matrix for efficiency.
    H = sparse(H);

end
