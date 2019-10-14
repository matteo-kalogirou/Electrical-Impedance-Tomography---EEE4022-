function [ avP ] = getAvPower( X )
% Returns the average power on each channel in a measurement set X
% WHere X is [Nreading * 16 channels * 16 frames]

frame_ampl = range(X)./2;                   % amplitude estimate
                                            % (sinusoidal signals)
frame_power = (frame_ampl.^2)/2;            % average signal power
avP         = mean(frame_power, 3);         % average power through frames


end

