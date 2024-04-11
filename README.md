
# Rate-Adaptive Protograph-Based Raptor-Like LDPC Code for Continuous-Variable Quantum Key Distribution

This file provides information on how to use the rate-adaptive protograph-based Raptor-like LDPC codes in continuous-variable quantum key distribution.

## How to use the code?
The code definitions are inside the "PCM" directory. The files starting with "H" are the parity check matrices (PCM) used in decoders.
 
Example MATLAB scripts on how to use the code are provided in the "scripts" directory. The code definition files have the ".mat" extension, but they can be also used in other programming languages/environments.

## Setup
To run the example scripts, MATLAB needs to be installed along with the Communications Toolbox.

## Naming Conventions
- "H_AZCW.mat" is the PCM used to simulate the all-zero codeword(AZCW). It is used mainly in Monte-Carlo simulations to get the performance of the code.
- For the syndrome decoding used in CV-QKD, "H_AZCW.mat" is enough, no generator matrix is necessary.

## To-Do
- Implementation of syndrome decoder


## Acknowledgment
This work was funded by the German Federal Ministry of Education and Research (BMBF) under grant agreement 16KISQ056 (DE-QOR).