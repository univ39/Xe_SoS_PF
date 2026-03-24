function [DATA] = noise_data(DATA_unnoised,SNR,noise_power)
%NOISE_DATA add gaussian distributed noise to the image, with defined SNR
%   SNR=mean(data)^2/std(noise)^2
if ~exist('noise_power','var');noise_power=[];end

if isempty(noise_power)
    sig_power=mean(abs(DATA_unnoised(:)).^2);
    noise_power=sig_power./SNR;
end

    real_noise=sqrt(noise_power/2)*randn(size(DATA_unnoised));
    imag_noise=sqrt(noise_power/2)*randn(size(DATA_unnoised)); %% complex
    %% noise the data
    DATA=DATA_unnoised + real_noise + 1i*imag_noise;
    % figure;plot(sum(abs(DATA_unnoised),1))
    % hold on;plot(sum(abs(DATA),1))
end

