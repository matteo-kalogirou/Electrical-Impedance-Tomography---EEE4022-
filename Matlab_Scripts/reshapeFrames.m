function [ shapedFrames ] = reshapeFrames( Frames )
% Removes three readings from the input matrix to conform to the
% conventional 16 electode model's measurement arrangement
%
% Input: nx16 matrix representing the 16 differential voltages obtained
% from the 32-electrode system, (interleaved voltage and current electrodes)
%
% Output: nx13 matrix representing the 13 possible voltage measurements
% obtainable from a 16 electrode system who's active electrode are occluded
% from the differential votlage measurements

if(size(Frames, 2) ~= 16)
    display('Error in input shape. Must be Nx16.');
    return;
end

shapedFrames = zeros(size(Frames, 1), 13, 16);

for i=1:16
    l = i-1; if(l<1) l=16; end
    r = i+1; if(r>16) r = 1; end;
    scrap = sort([l, i, r]);
    
    temp = Frames(:,:,i);
    for j=1:3
        temp(:,scrap(j)) = [];
        scrap = scrap-1;
    end
    
    shapedFrames(:,:,i) = temp;

end


end

