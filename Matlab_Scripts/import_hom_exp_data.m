%
% Import all homogenous trial data from /Results
% Average the results
%

results_path_300 = './Results/300Hz/hom_test';  % double delimitted by ';'
results_path_1k = './Results/1kHz/hom_test';  % double delimitted by ';'

for k=1:5
    % Convert all commas in files to decimal points
    
    for i=0:15
        comma2point_overwrite([results_path_300 num2str(k) '/frame' num2str(i) '.txt']);
        comma2point_overwrite([results_path_1k num2str(k) '/frame' num2str(i) '.txt']);
    end

    % Read in the files
    for i = 1:16
       Frames_300(:,:,i) =  dlmread([results_path_300 num2str(k) '/frame' num2str(i-1) '.txt'], ';', 0,0);
       Frames_1k(:,:,i) =  dlmread([results_path_1k num2str(k) '/frame' num2str(i-1) '.txt'], ';', 0,0);
    end

    % Calculate the RMS value of each channel
    for i=1:16
       RmsFrames_300(:,:,i) = rms(Frames_300(:,:,i), 1);
       RmsFrames_1k(:,:,i) = rms(Frames_1k(:,:,i), 1);
    end      
end

hom_exp_data_300 = reshape(RmsFrames_300, [256,1]);
hom_exp_data_1k = reshape(RmsFrames_1k, [256,1]);

% figure();
% plot(hom_data); axis tight;
% xlabel('Sample Number'); ylabel('RMS differential voltage');
% title('Homogenous Data');

clear RmsFrames_300 RmsFrames_1k Frames_300 Frames_1k results_path_300 results_path_1k i k

