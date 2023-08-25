clear all; close all; 
%% time based
% Nn+1=Nn/(1+dn)

initialLearnRate = 0.001; %0.001 si on initialise le réseau, sinon mettre un learningRate faible.
decayRate = 0.00001;

iteration=1:10:100000;


%learningRate = (initialLearnRate ./ (1+decayRate*iteration)).*((sin(iteration/64)+(1.1)));
learningRate = (initialLearnRate ./ (1+decayRate*iteration));
figure(111)
plot(iteration,learningRate)
title('time based')

%% step based

initialLearnRate = 0.1; %0.001 si on initialise le réseau, sinon mettre un learningRate faible.
decayRate = 0.95;
r=10;
iteration=1:1:1000;


learningRate = initialLearnRate*decayRate.^((1+iteration)./r);
figure(2)
plot(iteration,learningRate)
title('step based')

%% exponential

initialLearnRate = 0.1; %0.001 si on initialise le réseau, sinon mettre un learningRate faible.
decayRate = 0.005;

iteration=1:1:1000;

learningRate = initialLearnRate*exp(-decayRate.*iteration);
figure(3)

plot(iteration,learningRate)
title('exponential')

