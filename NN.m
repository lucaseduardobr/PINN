function parameters = NN(parameters,numLayers,numNeurons,numInput,numOutput,nameCouche)

%% Define Deep Learning Model

% Specify the number of layers and the number of neurons for each layer.



% Initialize the parameters for the first fully connect operation. The first fully connect operation has two input channels.

%parameters = struct;

sz = [numNeurons numInput]; %1
parameters.(nameCouche).("fc1").Weights = initializeHe(sz,numInput); % 1 ...
parameters.(nameCouche).("fc1").Bias = initializeZeros([numNeurons 1]);

%parameters.(nameCouche).("fc1").Bias = initializeC([numNeurons 1],0);
% Initialize the parameters for each of the remaining intermediate fully connect operations.

for layerNumber=2:numLayers-1
    name = "fc"+layerNumber;

    sz = [numNeurons numNeurons];
    numIn = numNeurons;
    parameters.(nameCouche).(name).Weights = initializeHe(sz,numIn);
    parameters.(nameCouche).(name).Bias = initializeZeros([numNeurons 1]);
    parameters.(nameCouche).(name).Alpha = initializeVariable([1 1],0.05);
    %parameters.(name).Bias = initializeC([numNeurons 1],0);
end

% Initialize the parameters for the final fully connect operation. The final fully connect operation has /one/ "two" output"s" channel.

sz = [numOutput numNeurons];
numIn = numNeurons;
parameters.(nameCouche).("fc"+ numLayers).Weights = initializeHe(sz,numIn);
parameters.(nameCouche).("fc"+ numLayers).Bias = initializeVariable([numOutput 1],1); 
parameters.(nameCouche).("fc"+ numLayers).Alpha = initializeVariable([1 1],0.05);
%parameters.(nameCouche).("fc"+ numLayers).Bias = initializeC([numOutput 1],1);
%parameters.(nameCouche).("fc" + numLayers).Bias = dlarray();