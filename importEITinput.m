function [ X ] = importEITinput( frame_folder )
% Accesses a folder of 16 frames of data and returns a vector of 1x256
% points to represent the experimental differential voltages on the tank
% Input: path-to-file eg ~/Users/matteokalogirou/Documents/
% files msut be store in this folder and must be .txt and named 'frameX.txt'


X = [];

for i=0:15
    comma2point_overwrite([frame_folder 'frame' num2str(i) '.txt']);
end

% Read in the files
for i = 1:16
   Frames(:,:,i) =  dlmread([frame_folder 'frame' num2str(i-1) '.txt'], ';', 0,0);
end

% Calculate the RMS value of each channel
for i=1:16
   RmsFrames(:,:,i) = rms(Frames(:,:,i), 1);
end

X = reshape(RmsFrames, [256,1]);

clear Frames RmsFrames;   


end

