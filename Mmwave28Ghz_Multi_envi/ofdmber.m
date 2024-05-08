close all;
clear all;
load('testcir.mat')
rng('shuffle');
Nsub = 1024;
Ncp = 1023;
SNRdB = 0.1:0.1:30;
n = 1;
M = 2^n;
BER_ls = zeros(size(SNRdB));
BER_mmse = zeros(size(SNRdB));
SNR = zeros(size(SNRdB));
s=0;
iter=100;
%%
for i=1:100
for L =1:iter
h = allcir{L};
h = h;
mes_bits = randi([0,M-1],[Nsub,1]); % generating random message
lenc = length(h) + Nsub +Ncp-1;
chnoise = (randn(1,lenc) + 1j*randn(1,lenc));
chnoises = chnoise/sqrt( mean(abs(chnoise).^2));    % nise generation
mes_mod  = pskmod(mes_bits,M);
mes_ifft = ifft(mes_mod,1024);
mes_cp   = [mes_ifft(Nsub-Ncp+1:Nsub);mes_ifft];
mes_pow  = mean(abs(mes_cp).^2);
for k = 1 : length(SNR)
%mes_tx   = (sqrt(snr).*mes_cp);
SNR(k)   = 10^(SNRdB(k)/10);
mes_tx   = (sqrt(SNR(k)).*mes_cp);
mes=mes_tx;
mes_tx   = reshape(mes_tx,[1,length(mes_tx)]);
sig_pow  = mean(abs(mes_tx).^2);
%%


channel = conv(h,mes_tx);

mes_rx  = channel +chnoise;
noise_pow = mean(abs(chnoises).^2);

sig_power = mean(abs(channel).^2);
%%
mes_rxwocp= mes_rx(Ncp + 1 : end);
mes_rxfft = fft(mes_rxwocp,1024);

pl = fft(mes_ifft,1024)'*fft(mes_ifft,1024);
est_ls    = conj(mes_rxfft .* fft(mes_ifft,1024)');
est_mmse = conj(est_ls)./(mean(abs(est_ls).^2)+1/SNR(k));
%%
est_mesmod_ls = mes_rxfft./est_ls;
est_mesmod_mmse = mes_rxfft.*est_mmse;
est_mesdemod_ls = pskdemod(est_mesmod_ls,M);
est_mesdemod_mmse = pskdemod(est_mesmod_mmse,M);
s=s+sum(reshape(est_mesdemod_ls,[length(est_mesdemod_ls),1]) ~= mes_bits);
BER_ls(k) = BER_ls(k) +sum(reshape(est_mesdemod_ls,[length(est_mesdemod_ls),1]) ~= mes_bits);
BER_mmse(k) = BER_mmse(k) +sum(reshape(est_mesdemod_mmse,[length(est_mesdemod_mmse),1]) ~= mes_bits);
end
end
end
BER_ls =BER_ls/(100*iter*Nsub);
BER_mmse =BER_mmse/(100* iter*Nsub);
semilogy(SNRdB,BER_ls)
hold on;
semilogy(SNRdB,BER_mmse);
axis tigth ;grid on;





%%