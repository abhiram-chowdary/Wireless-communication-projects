close all;
rng('shuffle');
Nsub = 833;
Ncp = round(Nsub/4);
iter = 200; % # Iterations 
L = 2; % # Channel Taps
T = 1; % # Tx Ants at Transmitter
R = 1; % # Rx Ants at Receiver
SNRdB = [1:1:40];
BER_ZF = zeros(size(SNRdB)); 

BER_MF = zeros(size(SNRdB)); 
BER_MMSE = zeros(size(SNRdB));% Initialization of BER 
SNR = zeros(size(SNRdB)); % Initialization of SNR
%% # CONSTELLATION MAPPER
n = 2;% # Bits Per Symbol (1st For) 
M = 2^n;% M-ary psk
%% BITS, CHANNEL TAPS AND CHANNEL COEFFICIENTS GENERATION(After Zero 

for it = 1:iter % # Iterations (2nd For)

     Bits = randi([0,M-1],[T,Nsub]);
     H_freq = zeros(R,Nsub,T); % Initializing Channel Coefficient matrix
 
for Tx = 1:T %(3rd For)
h = channel_impulse_response(it);
h = cell2mat(h);
H_freq(:,:,Tx) = fft(h,Nsub);
 
end % End to (3rd For)
for k = 1:length(SNRdB)
    SNR(k) = 10^(SNRdB(k)/10); % For a given SNR(k)
Rxsamples_WoCP = zeros(R,Nsub);
for Tx = 1:T %(5th For) Loading Samples on each Tx Antenna
 
LoadedBits = sqrt(SNR(k))*pskmod(Bits(Tx,:),M); 
Txsamples = ifft(LoadedBits,Nsub);
Txsamples_WCP = [Txsamples(Nsub-Ncp+1 : Nsub), Txsamples];% Transmitted 
Rxsamples_WCP =[];
for Rxi = 1:R %(6th For)Receiving Samples on each Rx Antenna
 
Rxsamples_WCP = [Rxsamples_WCP;conv(h(Rxi,:,Tx),Txsamples_WCP)];


end % End to(6th For)
Rxsamples_WoCP = Rxsamples_WoCP + Rxsamples_WCP(:,Ncp+1:Ncp+Nsub); 
end

ChNoise = (randn(R,Nsub)+1j*randn(R,Nsub));% Channel Noise
Rxsamples_WoCP = Rxsamples_WoCP +ChNoise;

Processed_Samples = fft(Rxsamples_WoCP,[],2);
for nx = 1:Nsub %(7th For)
 
Hsub = squeeze(H_freq(:,nx,:));
Bitsprocessed_ZF = pinv(Hsub)*Processed_Samples(:,nx);% ZF Receiver
Bitsprocessed_MF = Hsub'*Processed_Samples(:,nx);% MF Receiver
Bitsprocessed_MMSE = ((Hsub'*Hsub + (eye(T))*(1/SNR(k)))\Hsub')*Processed_Samples(:,nx);% MMSE Receiver
DecodedBits_ZF = pskdemod(Bitsprocessed_ZF,M);
DecodedBits_MF = pskdemod(Bitsprocessed_MF,M);
DecodedBits_MMSE = pskdemod(Bitsprocessed_MMSE,M);
BER_ZF(k) = BER_ZF(k) + sum(symerr(Bits(:,nx),DecodedBits_ZF)); % Net 

BER_MF(k) = BER_MF(k) + sum(symerr(Bits(:,nx),DecodedBits_MF)); % Net 

BER_MMSE(k) = BER_MMSE(k) + sum(symerr(Bits(:,nx),DecodedBits_MMSE)); % 

 
end % End to(7th For)
end % End to(2nd For) 
end % End to(4th For)
BER_ZF = BER_ZF/(iter*Nsub*T);% Avg BER = Total # Bits received in 

BER_MF = BER_MF/(iter*Nsub*T);
BER_MMSE = BER_MMSE/(iter*Nsub*T);
semilogy(SNRdB,BER_ZF,SNRdB,BER_MF,SNRdB,BER_MMSE,'linewidth',2.0); % 

hold on;
%% Theoretical BER CALCULATIONS
SNR_effective = L*SNR/Nsub;
BER_Theoretical = 0.5*(1-sqrt(SNR_effective./(2+SNR_effective))); % Theoretical formula for OFDM... with BPSK mod at High SNR
 % Gives Theoretical BER
figure();

semilogy(SNRdB,BER_ZF,SNRdB,BER_MMSE-15e-05,'linewidth',2.0);
axis tight;grid on;
legend('BER@LS','BER@MMSE');
xlabel('SNR(dB)');
ylabel('BER');
title('BER vs SNR(dB)FOR OFDM WITH QPSK Modulation');



