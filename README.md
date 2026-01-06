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

## Experiment-2 : 
2a.)Multiple Input Single Output: 
```matlab
%% Parameters 
Tsim       = 10;                 % total simulation time [s] 
lambda     = 179.69e-6;          % mean inter-arrival time [s] 
Npkt       = round(Tsim/lambda); % number of packets (~55660) 
 
Mset       = [1 2 4 8 16 32 64 128]; % number of Tx antennas to test 
Mmod       = 64;                % 64-QAM 
bitsPerSym = log2(Mmod);        % = 6 
bitsPerPkt = bitsPerSym;        % 1 QAM symbol per "packet" 
 
BW         = 10e6;              % 10 MHz 
NF_dB      = 5;                 % noise figure 
T0         = 290;               % noise temperature [K] 
kB         = 1.38064852e-23;    % Boltzmann 
Ptx_dBm    = 40;                % 40 dBm 
Ptx_W      = 10^(Ptx_dBm/10)/1000; % 10 W 
 
% Thermal noise power (W) in BW 
N0_W       = kB*T0*BW;          % thermal noise only 
% Add receiver NF 
N_W        = N0_W*10^(NF_dB/10); 
 
cohTime    = 10e-3;             % channel coherence time [s] 
pktRate    = 1/lambda;          % packets per second 
pktsPerBlk = max(1, round(cohTime * pktRate)); % packets per channel realization 
 
%% Preallocate BER results 
ber_all = zeros(length(Mset),1); 
 
%% Loop over number of Tx antennas 
for im = 1:length(Mset) 
   Mtx = Mset(im); 
   
   totBits  = 0; 
   bitErrs  = 0; 
   
   pktIdx   = 1; 
   while pktIdx <= Npkt 
       
       % New Rayleigh channel for this coherence block (1 x Mtx) 
       h = (randn(1,Mtx) + 1j*randn(1,Mtx))/sqrt(2); 
       
       % MRT beamforming vector (Mtx x 1), normalized 
       w = h'/norm(h);   % conjugate transpose implicitly in complex division 
       
       % Number of packets in this block (handle last partial block) 
       blkLen = min(pktsPerBlk, Npkt - pktIdx + 1); 
       
       % Generate bits for the whole block 
       dataBits = randi([0 1], bitsPerPkt*blkLen, 1); 
       
       % Map bits to 64-QAM symbols (unit average power) 
       txSym = qammod(dataBits, Mmod, 'InputType','bit', ... 
                      'UnitAveragePower', true); 
       
       % Scale transmit signal: x = w * sqrt(Ptx) * symbol 
       x = sqrt(Ptx_W) * (w * txSym.');        % (Mtx x blkLen) 
       
       % Received signal: y = h * x + noise 
       y_clean = h * x;                        % (1 x blkLen) 
       
       % AWGN noise 
       noise = sqrt(N_W/2) * (randn(1,blkLen) + 1j*randn(1,blkLen)); 
       y = y_clean + noise; 
       
       % Demodulate 64-QAM 
       rxBits = qamdemod(y.', Mmod, 'OutputType','bit', ... 
                         'UnitAveragePower', true); 
       
       % Count bit errors 
       bitErrs = bitErrs + sum(dataBits ~= rxBits); 
       totBits = totBits + numel(dataBits); 
       
       pktIdx = pktIdx + blkLen; 
   end 
   
   ber_all(im) = bitErrs / totBits; 
   fprintf('M=%d, BER=%g\n', Mtx, ber_all(im)); 
end 
 
%% Plot BER vs number of Tx antennas 
figure; 
semilogy(Mset, ber_all, 'o-','LineWidth',1.5); 
grid on; 
xlabel('Number of Tx antennas (M)'); 
ylabel('Bit Error Rate (BER)'); 
title('BER vs Number of Tx Antennas with MRT, 64-QAM, Flat Rayleigh'); 
```
2b.)Single Input Multiple Output :
```matlab
%% SIMO system BER simulation — Multiple Rx antennas, MRC combining 
 
% Parameters 
Tsim       = 10;                 % simulation time [s] 
lambda     = 179.69e-6;          % mean inter-arrival time [s] 
Npkt       = round(Tsim/lambda); % number of packets (~55660) 
 
RxSet      = [1 2 4 8 16 32 64 128]; % number of Rx antennas to test 
Mmod       = 64;                % 64-QAM 
bitsPerSym = log2(Mmod);        % = 6 
bitsPerPkt = bitsPerSym;        % 1 QAM symbol per "packet" 
 
BW         = 10e6;              % 10 MHz 
NF_dB      = 5;                 % noise figure 
T0         = 290;               % K 
kB         = 1.38064852e-23;    % Boltzmann 
Ptx_dBm    = 40;                % 40 dBm 
Ptx_W      = 10^(Ptx_dBm/10)/1000; % 10 W 
 
N0_W       = kB*T0*BW; 
N_W        = N0_W*10^(NF_dB/10); 
 
cohTime    = 10e-3;             
pktRate    = 1/lambda;          
pktsPerBlk = max(1, round(cohTime * pktRate)); 
 
ber_all = zeros(length(RxSet),1); 
 
for ir = 1:length(RxSet) 
   Nr = RxSet(ir); 
   
   totBits  = 0; 
   bitErrs  = 0; 
   pktIdx   = 1; 
   while pktIdx <= Npkt 
       % For each coherence block 
       h = (randn(Nr,1) + 1j*randn(Nr,1))/sqrt(2); % channel vector (Nr x 1) 
       
       blkLen = min(pktsPerBlk, Npkt - pktIdx + 1); 
 
       dataBits = randi([0 1], bitsPerPkt*blkLen, 1); 
       txSym    = qammod(dataBits, Mmod, 'InputType','bit', ... 
                         'UnitAveragePower', true); 
 
       x = sqrt(Ptx_W) * txSym.'; % transmitted symbol (1 x blkLen) 
       
       % Received signals, per antenna: (Nr x blkLen) 
       y_clean = h * x;                  
 
       % AWGN noise 
       noise = sqrt(N_W/2) * (randn(Nr,blkLen) + 1j*randn(Nr,blkLen)); 
       y = y_clean + noise; 
 
       % MRC combining: weigh by h*, sum antennas, normalize 
       y_combined = sum(conj(h) .* y,1) ./ sum(abs(h).^2); 
 
       rxBits = qamdemod(y_combined.', Mmod, 'OutputType','bit', ... 
                         'UnitAveragePower', true); 
 
       bitErrs = bitErrs + sum(dataBits ~= rxBits); 
       totBits = totBits + numel(dataBits); 
 
       pktIdx = pktIdx + blkLen; 
   end 
 
   ber_all(ir) = bitErrs / totBits; 
   fprintf('Nr=%d, BER=%g\n', Nr, ber_all(ir)); 
end 
 
figure; 
semilogy(RxSet, ber_all, 'o-','LineWidth',1.5); 
grid on; 
xlabel('Number of Rx antennas (Nr)'); 
ylabel('Bit Error Rate (BER)'); 
title('SIMO: BER vs Number of Rx antennas (MRC combining, 64-QAM, Rayleigh)');
```
