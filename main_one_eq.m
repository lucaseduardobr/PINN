close all; clearvars;clear all; clc;
%% loading data

load('Data_Exp_NoMass_MultiFreq_complex_phasecorrection');


%%% giving name to the variables %%%
mode = 6;
x = Beam.x';
w_exp = real(Beam.w(:,mode)) + j*(imag(Beam.w(:,mode)));
figure
plot(x,real(w_exp))
hold on
plot(x,imag(w_exp))
X_OBS_tot = x;
U_OBS_tot = w_exp;



%%% taking a part of the vector %%%%

n= length(w_exp);
debut = round(0.40*n)+1; % taking this position, the applied force does not "exist";
fin = round(0.80*n);
pas = 1;
X_OBS = x(debut:pas:fin);
U_OBS = w_exp(debut:pas:fin);
% plot(X_OBS,real(U_OBS),'x','LineWidth',20)
% hold on
% plot(X_OBS,imag(U_OBS),'o','LineWidth',15)


%%%% Normalizing %%%%
% I m normalizing each part (real and imaginary) with different constants
% then I m gonna to desnormalize in order to get a right residual
%normalize is important to help NN to converge since it doesn't explore
%large scale of data
%%% OBSERVATION: I m just normalizing  the data interval

aux2=1;
U=real(U_OBS);
Umax=max(abs(U));
U=U./Umax;

V=imag(U_OBS);
Vmax=max(abs(V));
V=V./Vmax;

%adding freq to the data struct
data.freq=Beam.freq(mode);


%%%%%%% Smoothing using Fourrier
% xa=X_OBS_tot(15:end-10);
% va=imag(U_OBS_tot(15:end-10));
deb=21;
fin=20;
%%% real
[~,u]=smoothing(X_OBS_tot(deb:end-fin),real(U_OBS_tot(deb:end-fin)),300);
%%% imag
[X_OBS,v]=smoothing(X_OBS_tot(deb:end-fin),imag(U_OBS_tot(deb:end-fin)),300);

U_OBS=u+j*v;

plot(X_OBS,real(U_OBS),'x','LineWidth',20)
hold on
plot(X_OBS,imag(U_OBS),'o','LineWidth',15)

%% getting a part of the values of the vector for futher testing
%it creates errors
% porc=0.2;   %%porcentage of the data located for testing
% 
% num_elementos = round(porc * length(X_OBS));
% indices = randperm(length(X_OBS), num_elementos);
% 
% X_Test = X_OBS(indices);
% U_OBS_Test=U_OBS(indices);
% %the permuta operation is changing the signal of the imaginary part
% %below I just fixing it, inverting the imaginary part
% U_OBS_Test=real(U_OBS_Test)+imag(U_OBS_Test)*-1j;
% X_OBS(indices) = [];
% U_OBS(indices) = [];

%ATTENTION U_TEST is determined by the neural network in further steps!




X_OBS_tot=X_OBS_tot';
U_OBS_tot=U_OBS_tot';

U_OBS_tot = real(U_OBS_tot) - imag(U_OBS_tot)*j;




%%% Generating random points


numInternalCollocationPoints = 1024; %6; %10000;
numTestPoints = 64;
pointSet = sobolset(1); 
points = net(pointSet,numInternalCollocationPoints); % pesca numInternalCollocationPoints dal poitSet
pointSetTest = sobolset(1);
pointsTest = net(pointSetTest,numTestPoints);
maxOBS = max(X_OBS);
minOBS = min(X_OBS);
dataX = (maxOBS-minOBS)*points(:,1)+minOBS; 
dataXTest = (maxOBS-minOBS)*pointsTest(:,1)+minOBS;
dataX = sort(dataX);

%five is just  to see the points
figure
y=5*ones(size(dataX));
plot(dataX,y,'o')
 hold on
dataX2=dataX(1:32:end);
%dataX2(end+1)=dataX(end); % to guarantee points on the extremity
dataX=dataX2;
y=5*ones(size(dataX));
plot(dataX,y,'x')



%Create an array datastore containing the training data.

ds = arrayDatastore([dataX]);


%% Define Deep Learning Model

%%% NN number 1 (of 2)
% Specify the number of layers and the number of neurons for each layer.

numLayers = 9;%9
numNeurons = 20;%20
numInput = 1;
numOutput = 2;
% Initialize the parameters for the first fully connect operation. The first fully connect operation has two input channels.

parameters = struct;
parameters = NN(parameters,numLayers,numNeurons,numInput,numOutput,"W");
% 
% numLayers = 9;
% numNeurons = 20;
% numInput = 2;
% numOutput = 2;
% parameters = NN(parameters,numLayers,numNeurons,numInput,numOutput,"E");

parameters.eConstante = initializeVariable([1 1],1.3);
parameters.eConstante2 = initializeVariable([1 1],1.3);
%% Load pre-trained neural network with network_save.m
load("network_save.mat") %la variable parameters est remplacée par la valeurs de parameters contenu dans le fichier 



%% Specify Training Options
%Train the model for 1000 epochs ... 250 %%% 3000 epochs with a mini-batch size of 1000.

numEpochs = 500000;
lossStop=1e-15;
miniBatchSize = 1024;

%To train on a GPU if one is available, specify the execution environment "auto". Using a GPU requires Parallel Computing Toolbox™ and a supported GPU device. For information on supported devices, see GPU Support by Release (Parallel Computing Toolbox) (Parallel Computing Toolbox).

executionEnvironment = "cpu";

%Specify ADAM optimization options.
%higher means more exploratory
initialLearnRate = 0.001; %0.001 si on initialise le réseau, sinon mettre un learningRate faible.

decayRate = 0.00005;
%decayRate = 0.00008;

%%
%Accelerate the model loss function using the dlaccelerate function. To learn more, see Accelerate Custom Training Loop Functions.

accfun = dlaccelerate(@modelLoss);


%% Train Network

%Train the network using a custom training loop.

mbq = minibatchqueue(ds, ...
    MiniBatchSize=miniBatchSize, ...
    MiniBatchFormat="BC", ...
    OutputEnvironment=executionEnvironment);


% Convert the COMSOL's observations to dlarray. For the input data points, specify format with dimensions "CB" (channel, batch).

X_OBS = dlarray(X_OBS,"CB");
U_OBS = dlarray(U_OBS);

%namalization constants

Umax = dlarray(Umax,"CB");
Vmax = dlarray(Vmax,"CB");


% If training using a GPU, convert the initial and conditions to gpuArray.

if (executionEnvironment == "auto" && canUseGPU) || (executionEnvironment == "gpu")
    X_OBS = gpuArray(X_OBS);
    U_OBS = gpuArray(U_OBS);
end
%Initialize the parameters for the Adam solver.

averageGrad = [];
averageSqGrad = [];


%Initialize the training progress plot.

% figure
% subplot(2,1,1)
% CO = colororder;
% lineLoss = animatedline(Color=CO(2,:));
% ylim([0 inf])
% xlabel("Iteration")
% ylabel("Loss")
% set(gca, 'YScale', 'log')
% grid on
% hold on 
% subplot(2,1,2)
% %CO = colororder;
% lineLossE1 = animatedline(Color=CO(2,:));
% hold on
% %CO = colororder;
% lineLossE2 = animatedline(Color=CO(3,:));
% hold off
% hold on
% %CO = colororder;
% lineLossE3 = animatedline(Color=CO(4,:));
% hold off
% ylim([0 inf])
% xlabel("Iteration")
% ylabel("E")
% legend('E1','E2','E3')
% set(gca, 'YScale', 'log')
% grid on

loss_history = NaN(4,numEpochs);

start = tic;

iteration = 0;
numTotPoints = numInternalCollocationPoints;


E=double((data.E1)/(10^11)*ones(size(dataX,1),1));
E_tot=(data.E1)/(10^11)*ones(size(X_OBS_tot,2),1);
X_Epoch=NaN(numEpochs/10,1);
E_Test=NaN(numEpochs/10,1);
E_Test_imag=NaN(numEpochs/10,1);
U_error=NaN(numEpochs/10,1);
U_error_imag=NaN(numEpochs/10,1);
Uxxxx_error=NaN(numEpochs/10,1);
Vxxxx_error_imag=NaN(numEpochs/10,1);

% options = trainingOptions('adam', ... % otimizador
%         'MaxEpochs',10, ... % número máximo de épocas
%         'MiniBatchSize',miniBatchSize, ... % tamanho do mini lote
%         'Plots','training-progress', ... % exibir gráficos de progresso
%        'OutputFcn',@(info) stopAtThreshold(info,9e-2)); % função de saída personalizada para verificar a perda

%name of the gif generate
name='loss.gif';
for epoch = 1:numEpochs
    %reset(mbq);
    %shuffle(mbq);
    %epoch
    %while hasdata(mbq)
        

        
        iteration = iteration + 1;
        
        %XX = next(mbq);
        % X = XX(1,:); % tocheck
        X = dlarray(dataX,"BC"); %quand il y a pas la Parallel Computing Toolbox
        XTest = dlarray(dataXTest,"BC");
        % Evaluate the model loss and gradients using dlfeval and the
        % modelLoss function.
        [loss,lossTest,gradients,Uxxxx,Uxxx,Uxx,Ux,Vxxxx,Vxxx,Vxx,Vx,A,B] = dlfeval(accfun,parameters,X,XTest,X_OBS,U_OBS,data,Umax,Vmax,numInternalCollocationPoints,numTestPoints);
        % Update learning rate.
        %learningRate = (initialLearnRate ./ (1+decayRate*iteration)).*((sin(iteration/128)+1.1));
        learningRate = initialLearnRate / (1+decayRate*iteration);
        %learningRate= 0.001;
        % Update the network parameters using the adamupdate function.
        [parameters,averageGrad,averageSqGrad] = ...
            adamupdate(parameters,gradients,averageGrad, ...
            averageSqGrad,iteration,learningRate);
        
       

   
        

    %end
    
    % RAR Method 
%     if mod(epoch,500)==0
%     
%         numTestPoints = 32; %6; %10000;
%         pointSet = sobolset(1); % general valori random from the sobol sequence ; in dimensione 1
%         pointsTest = net(pointSet,numTestPoints); % pesca numInternalCollocationPoints dal poitSet
%         dataXTest = (maxOBS-minOBS)*pointsTest(:,1)+minOBS; 
%         X = dlarray(dataXTest,"BC");
% 
%         indice = dlfeval(@modelLossRAR,parameters,X,data,w_max,numInternalCollocationPoints,loss);
%         
%         dataX = [dataX;
%                 dataXTest(indice)];
%         ds = arrayDatastore([dataX]);
%         
%         mbq = minibatchqueue(ds, ...
%             MiniBatchSize=miniBatchSize, ...
%             MiniBatchFormat="BC", ...
%             OutputEnvironment=executionEnvironment);
% 
%         numTotPoints = numTotPoints + length(indice);
%         fprintf('Nous avons ajouté : %d points où il faut satisfaire l EDP, nombre total de points : %d \n',length(indice),numTotPoints); 
%     end

    loss_tot = double(gather(extractdata(loss(1))));
    loss_pinn_real = double(gather(extractdata(loss(2))));
    loss_obs_real = double(gather(extractdata(loss(3))));    
    loss_obs_imag = double(gather(extractdata(loss(4)))); 

    % loss_tot = double(gather(extractdata(loss(1))));
    % loss_pinn_real = double(gather(extractdata(loss(2))));
    % loss_pinn_imag = double(gather(extractdata(loss(3))));
    % loss_obs_real = double(gather(extractdata(loss(4))));    
    % loss_obs_imag = double(gather(extractdata(loss(5)))); 


   % loss_tot = double(gather(extractdata(loss(1))));

    loss_history(:,epoch) = loss;
%     subplot(2,1,1)
%     addpoints(lineLoss,iteration, loss_tot(1));
    
%     title("Epoch: " + epoch + ", Elapsed: " + string(D) + ", Loss: " + loss_tot(1))
%     hold on
    
    
    
     fprintf('Epoch : %d, Loss pinn real : %d,  Loss obs real : %d ,Loss obs imag : %d , Loss total : %d\n',epoch,loss_pinn_real,loss_obs_real,loss_obs_imag,loss_tot); 
     %fprintf('Epoch : %d, Loss pinn real : %d, Loss pinn imag : %d, Loss obs real : %d ,Loss obs imag : %d , Loss total : %d\n',epoch,loss_pinn_real,loss_pinn_imag,loss_obs_real,loss_obs_imag,loss_tot); 
    
    
    
    if mod(epoch,10)==0
        
        D = duration(0,0,toc(start),Format="hh:mm:ss");
        X_Epoch(epoch/10)=epoch;
        U_Test = extractdata(model_tanh(parameters,dlarray(dataX',"CB"),"W"));
        %E_Test = extractdata(model_sigmoid(parameters,dlarray(dataX',"CB"),"E"));
        %E_Test = extractdata(model_tanh(parameters,dlarray(dataX',"CB"),"E"));
        E_Test(epoch/10) = extractdata(parameters.eConstante);
        E_Test_imag(epoch/10) = extractdata(parameters.eConstante2);
        
        %%% real
        U_OBS_Test_interpol_real = interp1(X_OBS_tot,real(U_OBS_tot)/dl2double(Umax),dataX,'spline');
        U_error(epoch/10)=(norm(U_Test(1,:)' - U_OBS_Test_interpol_real) / norm(U_OBS_Test_interpol_real))*100;
        %%% imaginary
        U_OBS_Test_interpol_imag = interp1(X_OBS_tot,imag(U_OBS_tot)/dl2double(Vmax),dataX,'spline');
        U_error_imag(epoch/10)=(norm(U_Test(2,:)' - U_OBS_Test_interpol_imag) / norm(U_OBS_Test_interpol_imag))*100;
        

        %%%%%% Ploting imaginary part and real part [predicted and exact]
        
        %%%real
        figure(100)
        subplot(2,1,1)      
        plot(X_OBS_tot,real(U_OBS_tot)/dl2double(Umax),dataX,(U_Test(1,:)),"x")
        xlabel('x')
        ylabel('w(x)')
        legend('true', 'predicted')    
        title("REAL , Epoch is "+epoch+" Error is "+sprintf('%.4f', U_error(epoch/10))+"%")
        
        %%%imaginary
        subplot(2,1,2)
        plot(X_OBS_tot,imag(U_OBS_tot)/dl2double(Vmax),dataX,(U_Test(2,:)),"x")
        xlabel('x')
        ylabel('w(x)')
        legend('true', 'predicted')    
        title("Imaginary , Epoch is "+epoch+" Error is "+sprintf('%.4f', U_error_imag(epoch/10))+"%")

        %%%saving
        if epoch == 10
        writeGIF(gcf,name,0)
        else
        writeGIF(gcf,name,1)
        end
        

        %plot(X_OBS_tot,E_function(X_OBS_tot)/1e11,dataX,(E_Test)',"x") %X_OBS_tot(1:end-2),E_COMSOL_tot,  ((extractdata(parameters.eConstante)).^2)
        
      

%         xlabel('Epochs')
%         ylabel('E')
%         legend('true', 'predicted')  
%         title("E.mean = " +sprintf('%.2f', mean(E_Test)) + ", Duration = " + string(D))
%         ylim([-inf,1])

        %%%%%% Ploting the error of the displacement
        

        figure(200)
        %U_Test = extractdata(model_tanh(parameters,dlarray(X_Test',"CB"),"W"));
        %%% real
        subplot(2,1,1)
        plot(X_Epoch,U_error)
        %set(gca, 'YScale', 'log')
        title("Real U error = " + sprintf('%.4f', U_error(epoch/10)) +"%  Duration = " + string(D) )
        xlabel('Epochs')
        ylabel('Error(%)')
        %%% imaginary
        subplot(2,1,2)
        plot(X_Epoch,U_error_imag)
        set(gca, 'YScale', 'log')
        title("Imaginary U error = " + sprintf('%.4f', U_error_imag(epoch/10)) +"%  Duration = " + string(D) )
        xlabel('Epochs')
        ylabel('Error(%)')
       
        


% 
%         E_error(epoch/10)=(norm(E_Test' - E_function(dataX)/1e11) / norm(E_function(dataX)/1e11))*100;
% 
%         plot(X_Epoch,E_error)
%         title("E error = "+ sprintf('%.2f', E_error(epoch/10)) + "%  Duration = " + string(D) )
%         xlabel('Epochs')
%         ylabel('Error(%)')

        %%% Saving
        if epoch == 10
        writeGIF(gcf,'error.gif',0)
        else
        writeGIF(gcf,'error.gif',1)
        end


        %%%%%%% Gradients Real

        aux2=2;
        figure(300)
        clf %erase the current plot
        subplot(4,1,1)
        
        Ux_true = gradient(real(U_OBS_tot(1:aux2:end))/dl2double(Umax),X_OBS_tot(1:aux2:end));
        plot(X_OBS_tot(1:aux2:end),Ux_true)  
        hold on
        
        plot(dataX,dl2double(Ux),'x') 
        hold off
        xlim([0.3 0.80])
        title("Ux   Duration = " + string(D))
        xlabel("x")

        subplot(4,1,2)
        
        Uxx_true = gradient(gradient(real(U_OBS_tot(1:aux2:end))/dl2double(Umax),X_OBS_tot(1:aux2:end)),X_OBS_tot(1:aux2:end));
        plot(X_OBS_tot(1:aux2:end),Uxx_true)  
        hold on
        
        plot(dataX,dl2double(Uxx),'x') 
        hold off
        xlim([0.3 0.8])
        title("Uxx")
        xlabel("x")
        subplot(4,1,3)
        Uxxx_true = gradient(gradient(gradient(real(U_OBS_tot(1:aux2:end))/dl2double(Umax),X_OBS_tot(1:aux2:end)),X_OBS_tot(1:aux2:end)),X_OBS_tot(1:aux2:end));
        plot(X_OBS_tot(1:aux2:end),Uxxx_true)  
        hold on
       
        plot(dataX,dl2double(Uxxx),'x') 
        hold off
        xlim([0.3 0.8])
        title("Uxxx")
        xlabel("x")

        subplot(4,1,4)

        %Uxxxx=dl2double(Uxxxx);
        Uxxxx_true = gradient(gradient(gradient(gradient(real(U_OBS_tot(1:aux2:end))/dl2double(Umax),X_OBS_tot(1:aux2:end)),X_OBS_tot(1:aux2:end)),X_OBS_tot(1:aux2:end)),X_OBS_tot(1:aux2:end));
        Uxxxx_interpol = interp1(X_OBS_tot(1:aux2:end),Uxxxx_true,dataX,'spline');



      

        try
        Uxxxx_error(epoch/10)=(norm(dl2double(Uxxxx)' - Uxxxx_interpol) / norm(Uxxxx_interpol))*100;
        catch
        Uxxxx_error(epoch/10)=(norm(dl2double(Uxxxx) - Uxxxx_interpol) / norm(Uxxxx_interpol))*100;    
        end
        
        plot(X_OBS_tot(1:aux2:end),Uxxxx_true)   
   
        hold on
        
        plot(dataX,dl2double(Uxxxx),'x')
        
        hold off
        xlim([0.3 0.8])
        title("Uxxxx Error= "+sprintf('%.2f', Uxxxx_error(epoch/10))+"%")
        xlabel("x")
        
        %%% saving

        if epoch == 10
        writeGIF(gcf,'gradient_real.gif',0)
        else
        writeGIF(gcf,'gradient_real.gif',1)
        end
        
        %%%%%%%%  Gradients Imaginary

        aux2=2;
        figure(310)
        clf %erase the current plot
        subplot(4,1,1)
        
        Vx_true = gradient(imag(U_OBS_tot(1:aux2:end))/dl2double(Vmax),X_OBS_tot(1:aux2:end));
        plot(X_OBS_tot(1:aux2:end),Vx_true)  
        hold on
        
        plot(dataX,dl2double(Vx),'x') 
        hold off
        xlim([0.3 0.8])
        title("Vx   Duration = " + string(D))
        xlabel("x")

        subplot(4,1,2)
        
        Vxx_true = gradient(gradient(imag(U_OBS_tot(1:aux2:end))/dl2double(Vmax),X_OBS_tot(1:aux2:end)),X_OBS_tot(1:aux2:end));
        plot(X_OBS_tot(1:aux2:end),Vxx_true)  
        hold on
        
        plot(dataX,dl2double(Vxx),'x') 
        hold off
        xlim([0.3 0.8])
        title("Vxx")
        xlabel("x")
        subplot(4,1,3)
        Vxxx_true = gradient(gradient(gradient(imag(U_OBS_tot(1:aux2:end))/dl2double(Vmax),X_OBS_tot(1:aux2:end)),X_OBS_tot(1:aux2:end)),X_OBS_tot(1:aux2:end));
        plot(X_OBS_tot(1:aux2:end),Vxxx_true)  
        hold on
       
        plot(dataX,dl2double(Vxxx),'x') 
        hold off
        xlim([0.3 0.8])
        title("Vxxx")
        xlabel("x")

        subplot(4,1,4)

        %Uxxxx=dl2double(Uxxxx);
        Vxxxx_true = gradient(gradient(gradient(gradient(imag(U_OBS_tot(1:aux2:end))/dl2double(Vmax),X_OBS_tot(1:aux2:end)),X_OBS_tot(1:aux2:end)),X_OBS_tot(1:aux2:end)),X_OBS_tot(1:aux2:end));
        Vxxxx_interpol = interp1(X_OBS_tot(1:aux2:end),Vxxxx_true,dataX,'spline');
        

        
        
        
        try
        Vxxxx_error_imag(epoch/10)=(norm(dl2double(Vxxxx)' - Vxxxx_interpol) / norm(Vxxxx_interpol))*100;
        catch
        Vxxxx_error_imag(epoch/10)=(norm(dl2double(Vxxxx)' - Vxxxx_interpol) / norm(Vxxxx_interpol))*100;    
        end   

        plot(X_OBS_tot(1:aux2:end),Vxxxx_true)   
        hold on 

        plot(dataX,dl2double(Vxxxx),'x')     
        hold off
        xlim([0.3 0.8])
        title("Vxxxx Imaginary Error= "+sprintf('%.2f', Vxxxx_error_imag(epoch/10))+"%")
        xlabel("x")






        %%% saving

        if epoch == 10
        writeGIF(gcf,'gradient_imag.gif',0)
        else
        writeGIF(gcf,'gradient_imag.gif',1)
        end
        

        %%%%%%%% Loss coeficients 
        
        figure(320)
        clf
        subplot(2,1,1)
        E_error=(norm(E_Test(epoch/10) - 0.7) / norm(0.7))*100;
        plot(X_Epoch,E_Test) 
        %hold on
        %plot(X_Epoch,ones(size(X_Epoch))*0.7)
        title("Coef Real: "+sprintf('%.4f', E_Test(epoch/10)) + " Epoch is "+epoch +" Duration = " + string(D) )
        xlabel('Epochs')
        ylabel('Coef real')

        subplot(2,1,2)
        E_error=(norm(E_Test_imag(epoch/10) - 0.7) / norm(0.7))*100;
        plot(X_Epoch,E_Test_imag) 
%       hold on
%       plot(X_Epoch,ones(size(X_Epoch))*0.7)
        title("Coef Imaginary: "+sprintf('%.4f', E_Test_imag(epoch/10)) + " Epoch is "+epoch +" Duration = " + string(D) )
        xlabel('Epochs') 
        ylabel('Coef imag')


        if epoch == 10
        writeGIF(gcf,'loss_coef.gif',0)
        else
        writeGIF(gcf,'loss_coef.gif',1)
        end


        %%%%%%%% Loss history 
        
        figure(330)
        clf
        x_loss = linspace(0,numEpochs,numEpochs);
        plot(x_loss,loss_history)
        legend("Total",'location','best')
        xlabel("epoch")
        ylabel("loss")
        set(gca, 'YScale', 'log')        
        legend('total loss', 'PDE loss real', 'PDE loss imag', 'w(x) real loss','w(x) imaginary loss')
        title("Loss = " + loss_history(1,end) + ", Durée = " + string(D))
                

        figure(400)
        clf
        Xpred=dl2double(X_OBS);
        A=dl2double(A);
        B=dl2double(B);
        plot(Xpred,A,'x')
        hold on
        plot(Xpred,B,'o')
        legend("NN prediction","Real value")

    end
    if loss_tot< lossStop
        disp("Loss was reached")
        break
    end
end

%% converting gif to mp4




%%
D = duration(0,0,toc(start),Format="hh:mm:ss");
drawnow
%Check the effectiveness of the accelerated function by checking the hit and occupancy rate.

accfun



%%%% Converting gif to mp4 files 
%read gif2mp4 to know what field means

gif2mp4('error')
gif2mp4('gradient_imag')
gif2mp4('gradient_real')
gif2mp4('loss','Loss_of_the_function')
gif2mp4('loss_coef','FrameRate',15)




%%Evaluate Model Accuracy 
%compare the predicted values of the deep learning model with the true solutions of the Poisson's equation using the l2 error.

% Save parameters into network_save.m file

save("network_save","parameters");

% Plot loss function
x_loss = linspace(0,numEpochs,numEpochs);
figure 
plot(x_loss,loss_history)
legend("Total")
xlabel("epoch")
ylabel("loss")
set(gca, 'YScale', 'log')
legend('total loss', 'PDE loss real', 'PDE loss imag', 'w(x) real loss','w(x) imaginary loss')
title("Loss = " + loss_history(1,end) + ", Durée = " + string(D))

%% Plotting prediction

%X_Test = X_OBS_tot(debut:fin); %(maxOBS-minOBS).*linspace(0,1,1000)+minOBS;
U_Test = extractdata(model_tanh(parameters,dlarray(X_Test',"CB"),"W"));
%E_Test = extractdata(model_sigmoid(parameters,dlarray(X_Test',"CB"),"E"));
%E_Test = extractdata(model_tanh(parameters,dlarray(X_Test',"CB"),"E"));


figure
%%% Real
subplot(2,1,1)
plot(X_OBS_tot,real(U_OBS_tot)/dl2double(Umax))
hold on 
err = norm(U_Test(1,:)' - real(U_OBS_Test)/dl2double(Umax)) / norm(real(U_OBS_Test)/dl2double(Umax));
plot(X_Test,real(U_OBS_Test)/dl2double(Umax),"x")
xlabel("x")
ylabel("w(x) real")
legend('true', 'predicted')
title("U real error "+ round(err*100,2) +" %" )

subplot(2,1,2)
plot(X_OBS_tot,imag(U_OBS_tot)/dl2double(Vmax))
hold on
plot(X_Test,imag(U_OBS_Test)/dl2double(Vmax),"x")
err = norm(U_Test(2,:)' - imag(U_OBS_Test)/dl2double(Vmax)) / norm(imag(U_OBS_Test)/dl2double(Vmax));
legend('true', 'predicted')
title("U imag error :  "+ round(err*100,2) +" %")
xlabel("x")
ylabel("w(x) imaginary")



%%

%save all the plots

saveas(100, 'displacement.png');
saveas(200, 'error.png');
saveas(300, 'real_gradient.png');
saveas(310, 'imag_gradient.png');
saveas(320, 'prediction.png');
saveas(330, 'loss.png');






% On regarde la différence sur la dérivée 4e 





 %Uxxxx_true = gradient(gradient(gradient(gradient(U_OBS_tot(1:aux:end),X_OBS_tot(1:aux:end)),X_OBS_tot(1:aux:end)),X_OBS_tot(1:aux:end)),X_OBS_tot(1:aux:end));



% Uxxxx_true= Uxxxx_true';


% Uxxxx_true = gradient(gradient(gradient(gradient(U_OBS_tot(1:aux:end),X_OBS_tot(1:aux:end)),X_OBS_tot(1:aux:end)),X_OBS_tot(1:aux:end)),X_OBS_tot(1:aux:end));
% 
% plot(X_OBS_tot(1:aux:end),Uxxxx_true);

%%Compare the gradients
%%% real 

figure(6)

try
scatter(extractdata(X),extractdata(Uxxxx),"x")
catch
scatter(extractdata(X),Uxxxx,"x")
end
hold on 
aux=32;

try
Uxxxx_true = deriveeOrdre4(real(U_OBS_tot(1:aux:end))'/dl2double(Umax),X_OBS_tot(1:aux:end));
catch
Uxxxx_true = deriveeOrdre4(real(U_OBS_tot(1:aux:end))/dl2double(Umax),X_OBS_tot(1:aux:end));  
end

plot(X_OBS_tot(1:aux:end),Uxxxx_true)
hold on
U = dl2double(model_tanh(parameters,dlarray(dataX',"CB"),"W"));
try
 Uxxxx_gradient=deriveeOrdre4(U(1,:),dataX);
catch
 Uxxxx_gradient=deriveeOrdre4(U(1,:)',dataX);
end

plot(dataX,Uxxxx_gradient)
hold off
%ylim([-1e8 1e8])
%xlim([0.3 0.55])
legend("Predicted","Exact","Ugradient")
xlim([0.3 0.8])
title("Real gradient")



%%% imaginary

figure(7)
try
scatter(extractdata(X),extractdata(Vxxxx),"x")
catch
scatter(extractdata(X),Vxxxx,"x")
end
hold on 

aux=32;

try
Vxxxx_true = deriveeOrdre4(imag(U_OBS_tot(1:aux:end))'/dl2double(Vmax),X_OBS_tot(1:aux:end));

catch
Vxxxx_true = deriveeOrdre4(imag(U_OBS_tot(1:aux:end))/dl2double(Vmax),X_OBS_tot(1:aux:end));  

end

plot(X_OBS_tot(1:aux:end),Vxxxx_true)


hold on

U = dl2double(model_tanh(parameters,dlarray(dataX',"CB"),"W"));


try
 Vxxxx_gradient=deriveeOrdre4(U(2,:),dataX);
catch
 Vxxxx_gradient=deriveeOrdre4(U(2,:)',dataX);

end


plot(dataX,Vxxxx_gradient)


hold off
%ylim([-1e8 1e8])
%xlim([0.3 0.55])
legend("Predicted","Exact","Ugradient")
xlim([0.3 0.8])
title("Imaginary")



%%% comparing displacements
%%%% real
figure(8)
subplot(2,1,1)
plot(dl2double(X_OBS(1:aux:end)),dl2double(real(U_OBS(1:aux:end))/dl2double(Umax)),'o')
hold on
plot(dataX,U(1,:),'x')
xlim([0.3 0.8])
legend("W_exact","W_pred")

%%%% imaginary
subplot(2,1,2)
plot(dl2double(X_OBS(1:aux:end)),dl2double(imag(U_OBS(1:aux:end))/dl2double(Vmax)),'o')
hold on
plot(dataX,U(2,:),'x')
xlim([0.3 0.8])
legend("W_exact","W_pred")


%% On clear l'accelerate function
clearCache(accfun)
% 
% function [stop,options] = stopAtThreshold(info, threshold)
%     % Verifica se a perda atingiu o limiar
%     stop = false;
%     if info.TrainingLoss(end) < threshold
%         stop = true;
%     end
%     
%     % Atualiza as opções de treinamento para interromper o treinamento
%     % quando a perda atingir o limiar
%     if stop
%         fprintf('A perda atingiu o limiar de %d O treinamento será interrompido',threshold );
%         options.MaxEpochs = info.Epoch;
%     end
% end



%err = norm(extractdata(U_Test) - U_OBS_tot) / norm(U_OBS_tot);

% subplot(2,1,1)
% plot(X_OBS_tot,U_OBS_tot)
% hold on 
% plot(extractdata(X_OBS),extractdata(U_OBS),"x")
% xlabel("x")
% ylabel("w(x)")
% legend('true', 'predicted')
% 
% subplot(2,1,2)
% plot(X_OBS_tot,2.1*1e11*ones(size(X_OBS_tot)))
% hold on 
% plot(extractdata(X_OBS),E_Test,"x")
% xlabel("x")
% ylabel("w(x)")