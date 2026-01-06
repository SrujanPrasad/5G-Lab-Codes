# 5G-Lab-Codes

## Experiment-1 : 
1a.) Generate sine wave : 
```matlab
clc;clear all;close all;
start_time=0;
end_time=1;
f=5;
fs=10000;
t=linspace(start_time,end_time,fs);
x=sin(2*pi*f*t);
plot(t,x);
xlabel('Time(s)');
ylabel('Amplitude');
title('Sine Wave');
grid on;
```

1c.)5G NR Waveform Generation : 
```matlab
% 5G NR downlink TX/RX + OFDM demod (simulation only, no USRP)
clear; clc;
%% 1) Configure 5G NR downlink carrier
cfg = nrDLCarrierConfig;      % Default downlink config (FR1)
% Generate time-domain NR waveform
[txWave, info] = nrWaveformGenerator(cfg);
% Get sample rate from waveform info
fs = info.ResourceGrids(1).Info.SampleRate;
%% 2) Simple channel (AWGN only)
snr = 30;                     % dB
rxSig = awgn(txWave, snr, 'measured');
%% 3) OFDM demodulation (correct use of nrOFDMDemodulate)
% Extract SCS carrier configuration
scsCarrier   = cfg.SCSCarriers{1};           % First SCS carrier
nrb          = scsCarrier.NSizeGrid;         % Number of resource blocks
scs          = scsCarrier.SubcarrierSpacing; % Subcarrier spacing (kHz)
initialNSlot = 0;                            % Start at slot 0
% Perform OFDM demodulation
% Syntax: rxGrid = nrOFDMDemodulate(waveform, nrb, scs, initialNSlot, ...)
rxGrid = nrOFDMDemodulate(rxSig, nrb, scs, initialNSlot, ...
                          'SampleRate', fs);
%% 4) Display info
disp('Size of rxGrid (subcarriers x OFDM symbols x Rx antennas):');
disp(size(rxGrid));
%% 5) Plot constellation of all REs (sanity check)
figure;
plot(real(rxGrid(:)), imag(rxGrid(:)), '.');
grid on;
xlabel('In-Phase');
ylabel('Quadrature');
title('Constellation of all received REs (no equalization/decoding)');
```
