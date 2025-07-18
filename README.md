# Channel Estimation and Equalization for PUSCH

## Overview
This project focuses on the channel estimation and equalization for the Physical Uplink Shared Channel (PUSCH) in 5G NR systems. The main goal is to improve the reliability of the communication by applying algorithms like Least Square (LS), Zero Force (ZF), and Minimum Mean Square Error (MMSE). These algorithms help in estimating the channel and performing equalization to combat channel impairments such as fading and interference.

## PUSCH Grid Use Case
The PUSCH is crucial for transmitting user data in 5G networks, particularly in the uplink direction. The grid of PUSCH resource elements is designed to allocate frequency and time resources dynamically based on traffic demands. This use case simulates the allocation and estimation of channels in this uplink shared resource.

## System Model
The system model includes:
- **Transmitter Model**: Converts binary data into modulated symbols and applies appropriate pilots for channel estimation.
- **Channel Model**: A wireless channel modeled by a tapped-delay line (TDL) system, which includes multi-path fading and Doppler shifts.
- **Receiver Model**: Receives the modulated signals, removes cyclic prefixes, performs FFT to convert time-domain signals to the frequency domain, estimates the channel, and equalizes the received signal.
![System Model](figures/system_model.png)
## Channel Estimation and Equalization Algorithms
The following algorithms are implemented to improve the performance of the PUSCH:
- **Least Square (LS)**: Direct estimation minimizing the error between the received and transmitted signals.
- **MMSE (Minimum Mean Square Error)**: A more sophisticated method that uses statistical information to reduce error in the estimation.
- **ZF (Zero Forcing)**: A simple equalization algorithm that inverts the estimated channel to recover the transmitted symbols, though it may amplify noise in low SNR conditions.

## Results
Simulation results show the effectiveness of channel estimation and equalization techniques in improving system performance. Key results include:
- **BER (Bit Error Rate)** improvement as the SNR (Signal to Noise Ratio) increases.
- **MSE (Mean Square Error)** reduction for the LS, ZF, and MMSE methods.
- A comparison of different modulation schemes (QPSK, 16-QAM, 64-QAM) demonstrates that MMSE outperforms other methods in terms of both BER and MSE.

Detailed charts and performance evaluations can be found in the project files.

## Installation
To get started, clone this repository:
```bash
git clone https://github.com/yourusername/Channel-Estimation-PUSCH.git
cd code
