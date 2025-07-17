# 📡 Mini Project Report: Channel Estimation and Equalization for PUSCH in 5G NR

This repository contains the mini project report titled:

> **"Designing Channel Estimation and Equalization Schemes for PUSCH using LS/ZF and MMSE Algorithms"**

---

## 📌 Project Overview

- **Program**: Viettel Digital Talent 2025 – Phase 1  
- **Student**: Nguyễn Khắc Kiên  
- **Supervisor**: Mr. Lương Xuân Hào  
- **Affiliation**: Faculty of Electronics and Telecommunications, University of Engineering and Technology, VNU

---

## 📖 Abstract

This project presents the design, simulation, and performance evaluation of channel estimation and equalization techniques for the Physical Uplink Shared Channel (PUSCH) in 5G New Radio (NR) systems. The focus is placed on implementing and comparing the effectiveness of LS (Least Squares), ZF (Zero Forcing), and MMSE (Minimum Mean Square Error) methods within a 5G uplink simulation framework.

---

## 📚 Key Contents

### 1. **Overview of 5G NR Architecture**
- 5G network components: gNB, CU/DU split, 5GC core
- RAN protocol stacks and transmission channel classification

### 2. **PUSCH Physical Layer Processing**
- Uplink signal chain: modulation, pilot insertion, IFFT, cyclic prefix
- Types of PUSCH resource allocation: dynamic vs configured grants

### 3. **System Model**
- Transmitter, TDL channel model, and receiver chain
- Resource grid modeling (OFDM slots, subcarriers)

### 4. **Channel Estimation Algorithms**
- **LS Estimator**: Simple and low-complexity; suffers under low SNR
- **MMSE Estimator**: Statistically optimal; higher accuracy, higher complexity

### 5. **Equalization Techniques**
- **Zero Forcing (ZF)**: Inverts channel response; amplifies noise
- **MMSE Equalizer**: Balances between noise suppression and inversion

---

## 📊 Simulation & Results

- **Metrics Evaluated**:
  - **BER** (Bit Error Rate)
  - **MSE** (Mean Square Error)
- **Scenarios Tested**:
  - Varying SNRs from 0 dB to 30 dB
  - Three modulation schemes: QPSK, 16-QAM, 64-QAM
- **Findings**:
  - MMSE consistently outperforms LS + ZF in both BER and MSE
  - Higher modulation schemes require more accurate estimation

---

## 📂 Files Included

- [`VDT_5G_Channel_Estimation.pdf`](./VDT_5G_Channel_Estimation.pdf) — Full report with figures, formulas, and data tables
- [`code/`](../code/) (optional) — Folder for MATLAB or Python simulation scripts (to be added)

---

## 🔍 Future Work

- Explore machine learning-based estimators for non-linear and time-varying 5G channels
- Implement real-time FPGA-based testing using SDR platforms
- Integrate SRS and advanced pilot design for enhanced estimation under mobility

---

## 📑 References

Key research and IEEE papers cited in the report, including recent works on deep learning channel estimation in OFDM systems. Full citation list in the PDF.

---

Thank you for reading!  
If you find this project useful, feel free to ⭐️ the repository or fork for your own research.
