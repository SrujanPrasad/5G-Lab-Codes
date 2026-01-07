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

## Experiment-2 : MIMO Beamforming 
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

MIMO Beamforming (example)
```matlab
%% Beamforming for MIMO-OFDM Systems
% This example shows how to model a point-to-point MIMO-OFDM system with
% beamforming. The combination of multiple-input-multiple-output (MIMO) and
% orthogonal frequency division multiplexing (OFDM) techniques have been
% adopted in recent wireless standards, such as 802.11x families, to
% provide higher data rate. Because MIMO uses antenna arrays, beamforming
% can be adopted to improve the received signal to noise ratio (SNR) which
% in turn reduces the bit error rate (BER).
%
% This example requires Communications Toolbox(TM).

%   Copyright 2014-2022 The MathWorks, Inc.

%% Introduction
% The term MIMO is used to describe a system where multiple transmitters or
% multiple receivers are present. In practice the system can take many
% different forms, such as single-input-multiple-output (SIMO) or
% multiple-input-single-output (MISO) system. This example illustrates a
% downlink MISO system. An 8-element ULA is deployed at the base station as
% the transmitter while the mobile unit is the receiver with a single
% antenna.
%
% The rest of the system is configured as follows. The transmitter power is
% 9 watts and the transmit gain is -8 dB. The mobile receiver is stationary
% and located at 2750 meters away, and is 3 degrees off the transmitter's
% boresight. An interferer with a power of 1 watt and a gain of -20 dB is
% located at 9000 meters, 20 degrees off the transmitter's boresight. 

% Initialize system constants
rng(2014);
gc = helperGetDesignSpecsParameters();

% Tunable parameters
tp.txPower = 9;           % watt
tp.txGain = -8;           % dB
tp.mobileRange = 2750;    % m
tp.mobileAngle = 3;       % degrees
tp.interfPower = 1;       % watt
tp.interfGain = -20;      % dB
tp.interfRange = 9000;    % m
tp.interfAngle =   20;    % degrees
tp.numTXElements = 8;       
tp.steeringAngle = 0;     % degrees
tp.rxGain = 108.8320 - tp.txGain; % dB

numTx= tp.numTXElements;

%%
% The entire scene can be depicted in the figure below.

helperPlotMIMOEnvironment(gc, tp);

%% Signal Transmission 
% First, configure the system's transmitter. 

[encoder,scrambler,modulatorOFDM,steeringvec,transmitter,...
    radiator,pilots,numDataSymbols,frmSz] = helperMIMOTxSetup(gc,tp);

%%
% There are many components in the transmitter subsystem, such as the
% convolutional encoder, the scrambler, the QAM modulator, the OFDM
% modulator, and so on. The message is first converted to an information
% bit stream and then passed through source coding and modulation stages to
% prepare for the radiation.

txBits = randi([0, 1], frmSz,1);
coded = encoder(txBits);
bitsS = scrambler(coded);
tx = qammod(bitsS,gc.modMode,'InputType','bit','UnitAveragePower',true);

%%
% In an OFDM system, the data is carried by multiple sub-carriers that are
% orthogonal to each other. 

ofdm1 = reshape(tx, gc.numCarriers,numDataSymbols);

%%
% Then, the data stream is duplicated to all radiating elements in the
% transmitting array

ofdmData = repmat(ofdm1,[1, 1, numTx]);
txOFDM = modulatorOFDM(ofdmData, pilots);
%scale
txOFDM = txOFDM * ...
    (gc.FFTLength/sqrt(gc.FFTLength-sum(gc.NumGuardBandCarriers)-1));

% Amplify to achieve peak TX power for each channel
for n = 1:numTx
    txOFDM(:,n) = transmitter(txOFDM(:,n));
end

%%
% In a MIMO system, it is also possible to separate multiple users spatial
% division multiplexing (SDMA). In these situations, the data stream is
% often modulated by a weight corresponding to the desired direction so
% that once radiated, the signal is maximized in that direction. Because in
% a MIMO channel, the signal radiated from different elements in an array
% may go through different propagation environments, the signal radiated
% from each antenna should be propagated individually. This can be achieved
% by setting CombineRadiatedSignals to false on the phased.Radiator
% component.

radiator.CombineRadiatedSignals = false;

%%
% To achieve precoding, the data stream radiated from each antenna in the
% array is modulated by a phase shift corresponding to its radiating
% direction. The goal of this precoding is to ensure these data streams add
% in phase if the array is steered toward that direction. Precoding can be
% specified as weights used at the radiator.

wR = steeringvec(gc.fc,[-tp.mobileAngle;0]);

%% 
% Meanwhile, the array is also steered toward a given steering angle, so
% the total weights are a combination of both precoding and the steering
% weights.

wT = steeringvec(gc.fc,[tp.steeringAngle;0]);
weight = wT.* wR;

%%
% The transmitted signal is thus given by

txOFDM = radiator(txOFDM,repmat([tp.mobileAngle;0],1,numTx),conj(weight));

%%
% Note that the transmitted signal, txOFDM, is a matrix whose columns
% represent data streams radiated from the corresponding elements in the
% transmit array.


%% Signal Propagation 
% Next, the signal propagates through a MIMO channel. In general, there are
% two propagation effects on the received signal strength that are of
% interest: one of them is the spreading loss due to the propagation
% distance, often termed as the free space path loss; and the other is the
% fading due to multipath. This example models both effects.

[channel,interferenceTransmitter,toRxAng,spLoss] = ...
    helperMIMOEnvSetup(gc,tp);
[sigFade, chPathG] =  channel(txOFDM);
sigLoss = sigFade/sqrt(db2pow(spLoss(1)));

%%
% To simulate a more realistic mobile environment, the next section also
% inserts an interference source. Note that in a wireless communication
% system, the interference is often a different mobile user.

% Generate interference and apply gain and propagation loss
numBits = size(sigFade,1);
interfSymbols = wgn(numBits,1,1,'linear','complex');
interfSymbols = interferenceTransmitter(interfSymbols);
interfLoss = interfSymbols/sqrt(db2pow(spLoss(2)));

%% Signal Reception
% The receiving antenna collects both the propagated signal as well as the
% interference and passes them to the receiver to recover the original
% information embedded in the signal. Just like the transmit end of the
% system, the receiver used in a MIMO-OFDM system also contains many
% stages, including OFDM demodulator, QAM demodulator, descrambler,
% equalizer, and Viterbi decoder.

[collector,receiver,demodulatorOFDM,descrambler,decoder] = ...
    helperMIMORxSetup(gc,tp,numDataSymbols);

rxSig = collector([sigLoss interfLoss],toRxAng);

% Front-end amplifier gain and thermal noise
rxSig = receiver(rxSig);

rxOFDM = rxSig * ...
    (sqrt(gc.FFTLength-sum(gc.NumGuardBandCarriers)-1)) / (gc.FFTLength);

% OFDM Demodulation
rxOFDM = demodulatorOFDM(rxOFDM);

% Channel estimation
hD = helperIdealChannelEstimation(gc,  numDataSymbols, chPathG);

% Equalization
rxEq = helperEqualizer(rxOFDM, hD, numTx);

% Collapse OFDM matrix
rxSymbs = rxEq(:);

rxBitsS = qamdemod(rxSymbs,gc.modMode,'UnitAveragePower',true,...
    'OutputType','bit');
rxCoded = descrambler(rxBitsS);
rxDeCoded = decoder(rxCoded);
rxBits = rxDeCoded(1:frmSz);

%% 
% A comparison of the decoded output with the original message stream
% suggests that the resulting BER is too high for a communication system.
% The constellation diagram is also shown below:

ber = comm.ErrorRate;
measures = ber(txBits, rxBits);
fprintf('BER = %.2f%%; No. of Bits = %d; No. of errors = %d\n', ...
    measures(1)*100,measures(3), measures(2));

%%

constdiag = comm.ConstellationDiagram('SamplesPerSymbol', 1,...
    'ReferenceConstellation', [], 'ColorFading',true,...
    'Position', gc.constPlotPosition);
% Display received constellation
constdiag(rxSymbs);

%%
% The high BER is mainly due to the mobile being off the steering direction
% of the base station array. If the mobile is aligned with the steering
% direction, the BER is greatly improved.

tp.steeringAngle = tp.mobileAngle;

% Steer the transmitter main lobe
wT = steeringvec(gc.fc,[tp.steeringAngle;0]);

[txBits, rxBits,rxSymbs] = helperRerunMIMOBeamformingExample(gc,tp,wT);

reset(ber);
measures = ber(txBits, rxBits);
fprintf('BER = %.2f%%; No. of Bits = %d; No. of errors = %d\n', ...
    measures(1)*100,measures(3), measures(2));

%%
constdiag(rxSymbs);

%% 
% Therefore, the system is very sensitive to the steering error. On the
% other hand, it is this kind of spatial sensitivity that makes SDMA 
% possible to distinguish multiple users in space.

%% Phase Shifter Quantization Effect
% The discussion so far assumes that the beam can be steered toward the
% exact desired direction. In reality, however, this is often not true,
% especially when the analog phase shifters are used. Analog phase shifters
% have only limited precision and are categorized by the number of bits
% used in phase shifts. For example, a 3-bit phase shifter can only
% represent 8 different angles within 360 degrees. Thus, if such
% quantization is included in the simulation, the system performance
% degrades, which can be observed from the constellation plot.

% analog phase shifter with quantization effect
release(steeringvec);
steeringvec.NumPhaseShifterBits = 4;
wTq = steeringvec(gc.fc,[tp.steeringAngle;0]);

[txBits, rxBits,rxSymbs] = helperRerunMIMOBeamformingExample(gc,tp,wTq);

reset(ber);
measures = ber(txBits, rxBits);
fprintf('BER = %.2f%%; No. of Bits = %d; No. of errors = %d\n', ...
    measures(1)*100,measures(3), measures(2));

constdiag = comm.ConstellationDiagram('SamplesPerSymbol', 1,...
    'ReferenceConstellation', [], 'ColorFading',true,...
    'Position', gc.constPlotPosition);
constdiag(rxSymbs);

%% Summary
% This example shows a system level simulation of a point-to-point
% MIMO-OFDM system employing beamforming. The simulation models many system
% components such as encoding, transmit beamforming, precoding, multipath
% fading, channel estimation, equalization, and decoding. 

%% Reference
% [1] Houman Zarrinkoub, Understanding LTE with MATLAB, Wiley, 2014
% 
% [2] Theodore S. Rappaport et al. Millimeter Wave Wireless Communications,
% Prentice Hall, 2014
```

## Experiment-3 : RR,MAX C/I AND PF
```matlab
clc; 
clear; 
 
% Network Parameters 
num_users = 10;               % Number of UEs 
BW = 10e6;                    % Total bandwidth [Hz] 
TTI = 1e-3;                   % Transmission Time Interval [s] 
num_TTI = 1000;               % Simulation duration 
PRBs = 100;                   % Number of PRBs 
PRB_BW = BW/PRBs;             % BW per PRB 
pkt_size = 1500*8;            % bits per packet (1500 bytes) 
 
% Channel Model (random SNR per user, changes every TTI) 
SNR_dB = 10 + 10*randn(num_users, num_TTI); % Average SNR fluctuates 
 
% User traffic: can be full-buffer or bursty (for this demo, use full-buffer) 
buffer = pkt_size * ones(num_users, 1); % Each user always ready to transmit 
 
% Results Storages 
throughput_RR      = zeros(num_users,1); 
throughput_MaxCI   = zeros(num_users,1); 
throughput_PF      = zeros(num_users,1); 
 
% --- ROUND ROBIN Scheduler --- 
for t = 1:num_TTI 
   active_users = 1:num_users;            % All UEs always active 
   user_idx = mod(t-1,num_users)+1;       % RR assignment 
   user_SNR = SNR_dB(:,t); 
   rate = PRB_BW * log2(1 + 10.^(user_SNR/10)); % per-user rate formula (bps) 
   throughput_RR(user_idx) = throughput_RR(user_idx) + rate(user_idx)*TTI; 
end 
 
% --- MAX C/I Scheduler --- 
for t = 1:num_TTI 
   user_SNR = SNR_dB(:,t); 
   [~,idx] = max(user_SNR);               % select user with max SNR 
   rate = PRB_BW * log2(1 + 10.^(user_SNR/10)); 
   throughput_MaxCI(idx) = throughput_MaxCI(idx) + rate(idx)*TTI; 
end 
 
% --- PROPORTIONAL FAIR Scheduler --- 
avg_rate_PF = zeros(num_users,1);          % running average user rate 
for t = 1:num_TTI 
   user_SNR = SNR_dB(:,t); 
   rate = PRB_BW * log2(1 + 10.^(user_SNR/10)); 
   pf_metric = rate ./ (avg_rate_PF + 1e-9); % PF metric 
   [~,idx] = max(pf_metric); 
   throughput_PF(idx) = throughput_PF(idx) + rate(idx)*TTI; 
   avg_rate_PF(idx) = 0.9*avg_rate_PF(idx) + 0.1*rate(idx); % update avg (slow window) 
end 
 
% Compute Jain's Fairness Index for each scheduler 
jain = @(x) (sum(x)^2) / (num_users * sum(x.^2)); 
 
fairness_RR = jain(throughput_RR); 
fairness_MaxCI = jain(throughput_MaxCI); 
fairness_PF = jain(throughput_PF); 
 
% Total System Throughput (Mbps) 
systemThr_RR = sum(throughput_RR) / (num_TTI*TTI)/1e6; 
systemThr_MaxCI = sum(throughput_MaxCI) / (num_TTI*TTI)/1e6; 
systemThr_PF = sum(throughput_PF) / (num_TTI*TTI)/1e6; 
 
% Display Results 
fprintf('--- Round Robin ---\nSystem Throughput: %.2f Mbps, Jain Fairness: %.3f\n', systemThr_RR, fairness_RR); 
fprintf('--- Max C/I ------\nSystem Throughput: %.2f Mbps, Jain Fairness: %.3f\n', systemThr_MaxCI, fairness_MaxCI); 
fprintf('--- Prop. Fair ---\nSystem Throughput: %.2f Mbps, Jain Fairness: %.3f\n', systemThr_PF, fairness_PF); 
 
% Plot Results 
figure; 
bar([throughput_RR throughput_MaxCI throughput_PF]/(num_TTI*TTI)/1e6); grid on; 
xlabel('User Index'); 
ylabel('Average Throughput per user (Mbps)'); 
legend('Round Robin','Max C/I','Prop. Fair'); 
title('User Throughput Comparison'); 
 
figure; 
bar([fairness_RR fairness_MaxCI fairness_PF]); 
set(gca,'xticklabel',{'RR','Max C/I','PF'}); 
ylabel('Jain''s Fairness Index'); 
title('Fairness Comparison'); 
```

## Experiment-4 : Path Loss Models : 
```matlab
function pathloss_expected_models()
    % Frequency in GHz
    fc = 3.5;
    % Distance in meters (log scale)
    d = logspace(0, 3.3, 200);  % 1 m to ~2000 m
    % Free-space reference path loss at 1 m
    PL0 = 32.4 + 20*log10(fc);
    % ---- Path-loss exponents ----
    % LOS
    n_RMa_LOS = 2.0;
    n_UMa_LOS = 2.7;
    n_UMi_LOS = 3.0;
    n_InH_LOS = 3.3;
    % NLOS
    n_RMa_NLOS = 2.7;
    n_UMa_NLOS = 3.2;
    n_UMi_NLOS = 3.5;
    n_InH_NLOS = 4.0;
    % ---- Compute LOS ----
    PL_RMa_LOS = PL0 + 10*n_RMa_LOS*log10(d);
    PL_UMa_LOS = PL0 + 10*n_UMa_LOS*log10(d);
    PL_UMi_LOS = PL0 + 10*n_UMi_LOS*log10(d);
    PL_InH_LOS = PL0 + 10*n_InH_LOS*log10(d);
    % ---- Compute NLOS ----
    PL_RMa_NLOS = PL0 + 10*n_RMa_NLOS*log10(d);
    PL_UMa_NLOS = PL0 + 10*n_UMa_NLOS*log10(d);
    PL_UMi_NLOS = PL0 + 10*n_UMi_NLOS*log10(d);
    PL_InH_NLOS = PL0 + 10*n_InH_NLOS*log10(d);
    % ---- Plot ----
    figure; hold on; grid on;
    % LOS (solid lines)
    semilogx(d, PL_RMa_LOS, 'g', 'LineWidth', 2);
    semilogx(d, PL_UMa_LOS, 'b', 'LineWidth', 2);
    semilogx(d, PL_UMi_LOS, 'r', 'LineWidth', 2);
    semilogx(d, PL_InH_LOS, 'm', 'LineWidth', 2);
    % NLOS (dashed lines)
    semilogx(d, PL_RMa_NLOS, 'g--', 'LineWidth', 2);
    semilogx(d, PL_UMa_NLOS, 'b--', 'LineWidth', 2);
    semilogx(d, PL_UMi_NLOS, 'r--', 'LineWidth', 2);
    semilogx(d, PL_InH_NLOS, 'm--', 'LineWidth', 2);
    xlabel('Distance (m)');
    ylabel('Path Loss (dB)');
    title('Expectation-Matched Path Loss Models (LOS + NLOS)');
    legend('RMa LOS','UMa LOS','UMi LOS','InH LOS', ...
           'RMa NLOS','UMa NLOS','UMi NLOS','InH NLOS', ...
           'Location','northwest');
    xlim([1 2000]);
    ylim([50 200]);
end
```

## Experiment-5 :Perfomance of OFDMA SU-MIMO
```matlab
% ============================================================
% OFDMA + SU-MIMO Throughput Simulation (SISO / 2x2 / 4x4)
% MATLAB version using qammod/qamdemod (no System Objects)
% ============================================================
clear; clc; close all;
%% ---------------- Simulation Parameters --------------------
numSubcarriers   = 256;         % OFDMA grid size
numSymbols       = 14;          % OFDM symbols per frame
bitsPerSymbol    = 4;           % 4 = 16QAM, 2 = QPSK
M                = 2^bitsPerSymbol;
SNRdB            = 0:2:30;      % SNR values to test
numFrames        = 150;         % Number of frames per SNR
antennaConfigs   = [1 1; 2 2; 4 4];   % SISO, 2x2, 4x4
%% For Gray-coded QAM using qammod/qamdemod
modOrder = M;
symbolMap = 0:modOrder-1;
%% ---------------- Result Storage ---------------------------
throughput = zeros(size(antennaConfigs,1), length(SNRdB));
% ============================================================
%                     MAIN SIMULATION LOOP
% ============================================================
for cfg = 1:size(antennaConfigs,1)
    Nt = antennaConfigs(cfg,1);  % Tx antennas
    Nr = antennaConfigs(cfg,2);  % Rx antennas
    fprintf("\nSimulating %dx%d SU-MIMO...\n", Nt, Nr);
    for si = 1:length(SNRdB)
        snr = SNRdB(si);
        noiseVar = Nt./(10^(snr/10));  % scaling for MIMO system
        bitsTotal = 0;
        bitsCorrect = 0;
        for f = 1:numFrames
            % ---------------------------------------------------------
            % 1. OFDMA Grid → bits → QAM symbols
            % ---------------------------------------------------------
            numBits = numSubcarriers * numSymbols * bitsPerSymbol * Nt;
            txBits  = randi([0 1], numBits, 1);
            % Bits → integers
            txInts = bi2de(reshape(txBits, bitsPerSymbol, []).','left-msb');
            % QAM modulate (average power normalized)
            txSymbols = qammod(txInts, M, 'gray', 'UnitAveragePower', true);
            % Reshape into OFDMA grid: SC × SYM × Nt
            txSymbols = reshape(txSymbols, numSubcarriers, numSymbols, Nt);
            % ---------------------------------------------------------
            % 2. Rayleigh MIMO Channel + AWGN
            % ---------------------------------------------------------
            H = (randn(Nr,Nt) + 1j*randn(Nr,Nt)) / sqrt(2);  % i.i.d Rayleigh
            rxSymbols = zeros(numSubcarriers, numSymbols, Nr);
            for sc = 1:numSubcarriers
                for sy = 1:numSymbols
                    txVec = squeeze(txSymbols(sc,sy,:));
                    noise = sqrt(noiseVar/2)*(randn(Nr,1)+1j*randn(Nr,1));
                    rxSymbols(sc,sy,:) = H * txVec + noise;
                end
            end
            % ---------------------------------------------------------
            % 3. MMSE Equalization
            % ---------------------------------------------------------
            W = (H'*H + noiseVar*eye(Nt)) \ H';   % Nt × Nr
            eqSymbols = zeros(numSubcarriers, numSymbols, Nt);
            for sc = 1:numSubcarriers
                for sy = 1:numSymbols
                    r = squeeze(rxSymbols(sc,sy,:));
                    eqSymbols(sc,sy,:) = W * r;
                end
            end
            % ---------------------------------------------------------
            % 4. QAM Demodulation
            % ---------------------------------------------------------
            eqVec = eqSymbols(:);
            rxInts = qamdemod(eqVec, M, 'gray', 'UnitAveragePower', true);
            % integers → bits
            rxBits = reshape(de2bi(rxInts, bitsPerSymbol, 'left-msb').', [], 1);
            % ---------------------------------------------------------
            % 5. Throughput Measurement
            % ---------------------------------------------------------
            bitsTotal   = bitsTotal + numBits;
            bitsCorrect = bitsCorrect + sum(rxBits == txBits);
        end % frames
        ber = 1 - bitsCorrect/bitsTotal;
        throughput(cfg,si) = bitsCorrect / numFrames; % bits/frame
        fprintf("  SNR = %2d dB → BER = %.3g\n", snr, ber);
    end % SNR
end % antenna configs
% ============================================================
%                     PLOTTING RESULTS
% ============================================================
figure; hold on; grid on;
plot(SNRdB, throughput(1,:), '-o', 'LineWidth', 2);
plot(SNRdB, throughput(2,:), '-s', 'LineWidth', 2);
plot(SNRdB, throughput(3,:), '-^', 'LineWidth', 2);
xlabel('SNR (dB)'); ylabel('Throughput (bits/frame)');
title('Throughput vs SNR for SISO / 2×2 / 4×4 SU-MIMO');
legend('SISO (1x1)', '2×2 SU-MIMO', '4×4 SU-MIMO', 'Location','northwest');
```
