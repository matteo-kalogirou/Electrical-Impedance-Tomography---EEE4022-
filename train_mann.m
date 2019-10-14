%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Script used to train a multiple neural network
%   Matteo Kalogirou
%   
%   Tomasz Rymarczyk 1,2 , Grzegorz K?osowski 3,* , Edward Koz?owski 3 and Pawe? Tch?rzewski 2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;

x = nn_input;
t = nn_output;

trainFcn = 'trainlm';
hiddenLayerSize = 20;

net = fitnet(hiddenLayerSize, trainFcn);

net.input.processFcns = {'removeconstantrows', 'mapminmax'};
net.output.processFcns = {'removeconstantrows', 'mapminmax'};

net.divideFcn = 'dividerand';
net.divideMode = 'sample';
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;

net.performFcn = 'mse';

N = 256;
%%
parfor i=1:N
    
   y = t(i, :);
   [output_nets{i}, ~] = train(net, x, y);
   
end


