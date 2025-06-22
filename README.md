# Digital Communication System Simulator in C++

This project is a C++ simulation of a complete digital communication system. It is designed to model and analyze the process of transmitting digital data over a noisy channel, implementing key stages such as error correction coding, modulation, demodulation, and decoding. The primary goal is to evaluate the system's performance by calculating the Bit Error Rate (BER) under various noise conditions.

## Key Features

- **Error Correction Coding:** Implements **Hamming(7,4)** code for forward error correction (FEC), capable of correcting single-bit errors.
- **Digital Modulation Schemes:** Includes three fundamental modulation techniques:
  - Amplitude-Shift Keying (ASK)
  - Phase-Shift Keying (PSK)
  - Frequency-Shift Keying (FSK)
- **Noisy Channel Simulation:** Models a communication channel with two types of noise:
  - Additive White Gaussian Noise (AWGN)
  - Multiplicative (fading-like) noise
- **Coherent Demodulation:** Performs demodulation by integrating the received signal with a reference carrier wave
- **Performance Analysis:** Calculates the **Bit Error Rate (BER)** by comparing the transmitted and received bitstreams to quantify the system's reliability

![[Pasted image 20250622040818.png]]
![[Pasted image 20250622040857.png]]
![[Pasted image 20250622040911.png]]
## Technologies & Libraries

- **C++17**
- **Eigen:** For linear algebra operations
- **matplotlib-cpp:** A C++ wrapper for Python's Matplotlib, used for plotting and visualization during development

## Configuration

The main simulation parameters can be configured directly in the `main()` function:
- `in_bits`: The input bit string to be transmitted.
- `alfa`: Controls the intensity of the additive noise.
- `beta`: Controls the intensity of the multiplicative noise.
- `config`: Sets the order in which noise is applied (1 or 2).
- You can also switch between `ASK`, `PSK`, and `FSK` modulation within the `transSystem` function
