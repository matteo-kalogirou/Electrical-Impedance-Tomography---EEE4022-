%
% Multiply Artificial Neural Network
% Trying to address the ill-defined problem by increasing the number of
% neural networks such that one network is assigned to each element_data of
% the output. Each netowrk still requires all 208 inputs but produces a
% single output
%
% Matteo Kalogirou
%


x = nn_input;       % 208 x num_trials
t = nn_output;      % 576 x num_trials

n_trials = size(nn_input, 2);
n_outputElements = size(nn_output, 1);

trainFcn = 'trainbr';
hiddenLayerSize = 100;

%%
net = cell(n_outputElements,1);
%Create 576 networks - one for each output element_data
for i=1:n_outputElements
    net{i} = fitnet(hiddenLayerSize, trainFcn);
    net{i}.divideFcn = 'dividerand';
    net{i}.divideMode = 'sample';
    net{i}.divideParam.trainRatio = 70/100;
    net{i}.divideParam.valRatio = 15/100;
    net{i}.divideParam.testRatio = 15/100;
    net{i}.performFcn = 'mse';
end

%% Train the networks
tnet = cell(n_outputElements,1);
for i=1:n_outputElements
   % apply all inputs to train for one element
   tnet{i} = train(net{i}, x, t(i,:));
end



