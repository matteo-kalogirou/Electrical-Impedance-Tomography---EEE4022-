
%
% Quantifying the noise in the system
%
%

clc;
clear;

%---MEAS 1
noise_path = './Results/0_input/hom_test_1/';

for i=0:15
    comma2point_overwrite([noise_path 'frame' num2str(i) '.txt']);
end
% Read in the files
for i = 1:16
   noise_1.Frames(:,:,i) =  dlmread([noise_path 'frame' num2str(i-1) '.txt'], ';', 0,0);
end
% Calculate the RMS value of each channel
for i=1:16
   noise_1.RmsFrames(:,:,i) = rms(noise_1.Frames(:,:,i), 1);
end

noise_1.rmsAmp = reshape(noise_1.RmsFrames, [256,1]);

%Get the mean and std of each channel for each frame
noise_1.FramesMean = mean(noise_1.Frames);
noise_1.meanMean = mean(noise_1.FramesMean, 3);
noise_1.FramesVar = var(noise_1.Frames);
noise_1.meanVar = mean(noise_1.FramesVar, 3);
noise_1.FramesSTD = std(noise_1.Frames);
noise_1.meanSTD = mean(noise_1.FramesSTD, 3);


%---MEAS 2
noise_path = './Results/0_input/hom_test_2/';

for i=0:15
    comma2point_overwrite([noise_path 'frame' num2str(i) '.txt']);
end
% Read in the files
for i = 1:16
   noise_2.Frames(:,:,i) =  dlmread([noise_path 'frame' num2str(i-1) '.txt'], ';', 0,0);
end
% Calculate the RMS value of each channel
for i=1:16
   noise_2.RmsFrames(:,:,i) = rms(noise_2.Frames(:,:,i), 1);
end

noise_2.rmsAmp = reshape(noise_2.RmsFrames, [256,1]);

%Get the mean and std of each channel for each frame
noise_2.FramesMean = mean(noise_2.Frames);
noise_2.meanMean = mean(noise_2.FramesMean, 3);
noise_2.FramesVar = var(noise_2.Frames);
noise_2.meanVar = mean(noise_2.FramesVar, 3);
noise_2.FramesSTD = std(noise_2.Frames);
noise_2.meanSTD = mean(noise_2.FramesSTD, 3);


% ---MEAS 3
noise_path = './Results/0_input/hom_test_3/';

for i=0:15
    comma2point_overwrite([noise_path 'frame' num2str(i) '.txt']);
end
% Read in the files
for i = 1:16
   noise_3.Frames(:,:,i) =  dlmread([noise_path 'frame' num2str(i-1) '.txt'], ';', 0,0);
end
% Calculate the RMS value of each channel
for i=1:16
   noise_3.RmsFrames(:,:,i) = rms(noise_3.Frames(:,:,i), 1);
end

noise_3.rmsAmp = reshape(noise_3.RmsFrames, [256,1]);

%Get the mean and std of each channel for each frame
noise_3.FramesMean = mean(noise_3.Frames);
noise_3.meanMean = mean(noise_3.FramesMean, 3);
noise_3.FramesVar = var(noise_3.Frames);
noise_3.meanVar = mean(noise_3.FramesVar, 3);
noise_3.FramesSTD = std(noise_3.Frames);
noise_3.meanSTD = mean(noise_3.FramesSTD, 3);

% ----SIGNAL
signal_path = './NewResults/hom_1k_1vpp_1/';

for i=0:15
    comma2point_overwrite([noise_path 'frame' num2str(i) '.txt']);
end
% Read in the files
for i = 1:16
   signal.Frames(:,:,i) =  dlmread([signal_path 'frame' num2str(i-1) '.txt'], ';', 0,0);
end
% Calculate the RMS value of each channel
for i=1:16
   signal.RmsFrames(:,:,i) = rms(signal.Frames(:,:,i), 1);
end
signal.rmsAmp = reshape(signal.RmsFrames, [256,1]);

%% Average the channels means and variances together


av_mean = mean([noise_1.meanMean; noise_2.meanMean; noise_3.meanMean], 1);
av_var = mean([noise_1.meanVar; noise_2.meanVar; noise_3.meanVar], 1);
av_std = mean([noise_1.meanSTD; noise_2.meanSTD; noise_3.meanSTD], 1);




%% Plot normal dist for each channel from one measurement frame

subject = noise_1.Frames(:,:,16);
x = -0.5:0.001:0.4;

figure(2);

for i =1:16   
    ch_pd{i} = fitdist(subject(:,i), 'Normal');        
    hold on;
    y = normpdf(x, ch_pd{i}.mu, ch_pd{i}.sigma);
    %Normalise
    y = y./max(y);
    plot(x, y);    
    text(ch_pd{i}.mu, y(find(y>=max(y)-0.06*i,1)), ['\leftarrow Ch_{' num2str(i) '}'], 'FontSize', 14);
    hold off;
end
set(gca, 'FontSize', 14);
xlabel('Noise signal voltage');


%% Normal distribution of means
x = -1:0.001:1;

pd = fitdist(av_mean', 'Normal');

plot(x, normpdf(x, pd.mu, pd.sigma));
% title('Normal Distribution of the mean noise signal in each channel.');
xlabel('Mean( mean (channel voltage))')
% xlabel(['$\overline{\overline{x}(channel voltage)}$'],'interpreter','latex','FontSize',14);
text(pd.mu, 2.35, '\downarrow \mu = -0.0047', 'FontSize', 14);
text(pd.sigma, 1.4, '\leftarrow \sigma = 0.1726', 'FontSize', 14);
set(gca, 'FontSize', 14);

%% Normal distribution of variance
x = -0.0001:0.000001:0.0001;

pd2 = fitdist(av_var', 'Normal');

plot(x, normpdf(x, pd2.mu, pd2.sigma));
% title('Normal Distribution of the average noise signal variance in each channel.');
% xlabel(['$\overline{\sigma(channel voltage)}$'],'interpreter','latex','FontSize',14);
xlabel('Mean( var( channel voltage))');
text(pd2.mu, 3.6e4, '\downarrow \mu = 7.87e-06', 'FontSize', 14);
text(pd2.sigma, 3.3e4, '\leftarrow \sigma = 1.13e-05', 'FontSize', 14);
axis([-1e-4 1e-4 0 4e4]);
set(gca, 'FontSize', 14);

%% Get non-zero signal measurement and calculate the SNR

figure(12);
bar(noise_1.meanMean);
xlabel('Output Channel');
ylabel('DC Offset [V]');
set(gca, 'FontSize', 14);
grid on;
axis([0 17 -0.5 0.5])


%% Frequency Analysis

Fs = 5000;
% x = signal.Frames(:,1,1);
figure();
f = [];
for j = 1:16
    for i =1:16
        x = noise_1.Frames(:,i,j);
        N = length(x);
        xdft = fft(x);
        xdft = xdft(1:N/2+1);
        psdx = (1/(Fs*N)) * abs(xdft).^2;
        psdx(2:end-1) = 2*psdx(2:end-1);
        freq = 0:Fs/length(x):Fs/2;    
        f = [f; freq];
    end        
end

f = sum(f,1);

plot(f, 10*log10(psdx));
grid on; axis tight;
set(gca, 'FontSize', 14);
xlabel('Frequency [Hz]');
ylabel('Power Spectral Desnity [dB/Hz]')


%% SNR for each channel - use RMS amplitude


snr1 = 20*log10( rms(reshape(signal.rmsAmp,[16,16]),1)./ rms(reshape(noise_1.rmsAmp,[16,16]),1));
% stem(snr1);
bar(snr1); grid on;
xlabel(' Output Channel');
ylabel('SNR [dB]');
set(gca, 'FontSize', 14);


%% Get the signal power at each channel then average through each frame then calculate SNR

noise_1.power = getAvPower(noise_1.Frames);
noise_2.power = getAvPower(noise_2.Frames);
noise_3.power = getAvPower(noise_3.Frames);
signal.power = getAvPower(signal.Frames);

snr = [ 10*log10(signal.power./noise_1.power);
        10*log10(signal.power./noise_2.power);
        10*log10(signal.power./noise_3.power); ];

snr = mean(snr); 
meanSNR = mean(snr);

bar(snr1); refline(0, meanSNR);
grid on; axis ([0 17 -1 30]);
xlabel(' Output Channel');
ylabel('SNR [dB]');
set(gca, 'FontSize', 14);



min






