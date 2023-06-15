function U = model_tanh(parameters,XX,nameCouche)


X = [XX];
numLayers =numel(fieldnames(parameters.(nameCouche)));

% first layer First fully connect operation.
weights = parameters.(nameCouche).fc1.Weights;
bias = parameters.(nameCouche).fc1.Bias;
U = fullyconnect(X,weights,bias);

% mid layers tanh and fully connect operations for remaining layers.
for i=2:numLayers-1
    name = "fc" + i;
    alpha = parameters.(nameCouche).(name).Alpha;
    U = tanh(20*alpha.*U);

    weights = parameters.(nameCouche).(name).Weights;
    bias = parameters.(nameCouche).(name).Bias;
    U = fullyconnect(U, weights, bias);
end

%last layer 
name = "fc" + numLayers;
alpha = parameters.(nameCouche).(name).Alpha;
U = tanh(20*alpha.*U);
weights = parameters.(nameCouche).(name).Weights;
bias = parameters.(nameCouche).(name).Bias;
U = fullyconnect(U, weights, bias);



end
