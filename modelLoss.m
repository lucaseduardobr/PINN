function [loss,lossTest,gradients,Uxxxx,Uxxx,Uxx,Ux,Vxxxx,Vxxx,Vxx,Vx] = modelLoss(parameters,X,XTest,X_OBS,U_OBS,data,Umax,Vmax,numInternalCollocationPoints,numTestPoints)%,C0

% Make predictions

U= model_tanh(parameters,X,"W");
%E = (model_tanh(parameters,X,"E"));
% E = (model_sigmoid(parameters,X,"E"));
U_pred_OBS = model_tanh(parameters,X_OBS,"W");
%ebis = (parameters.eConstante).*ones(size(X));
Er =(parameters.eConstante);
Ei=(parameters.eConstante2);
% Assign variables from data 
%e = ebis;
%L = data.L;
rho = data.rho;
b = data.b;
e = data.e;
S = b*e;
% Sbis = b*ebis;
%nu = data.nu;
I = b*e.^3/(12); 
% Ibis = b*ebis.^3/12;
%E = data.E1;
F = data.F;
freq = data.freq;
w = 2*pi*freq;
x_F_min = 0.1;
x_F_max = 0.2;




% UTest= model_tanh(parameters,XTest,"W");
% bpredTest = model_tanh(parameters,XTest,"b");
% 
% % On calcule le résidu avec des XTest, c'est à dire qu'on ne fera pas de
% % mise à jour à partir des résultats de ces points.
% gradientsUTest = dlgradient(sum(UTest,"all"),{XTest},EnableHigherDerivatives=true);
% UxTest = gradientsUTest{1};
% 
% gradientsUxTest = dlgradient(sum(UxTest,"all"),{XTest},EnableHigherDerivatives=true);
% UxxTest = gradientsUxTest{1};
% 
% gradientsUxxTest = dlgradient(sum(EpredTest.*UxxTest,"all"),{XTest},EnableHigherDerivatives=true);
% UxxxTest = gradientsUxxTest{1};
% 
% gradientsUxxxTest = dlgradient(sum(UxxxTest,"all"),XTest,EnableHigherDerivatives=true);
% UxxxxTest = gradientsUxxxTest;
% 
% %%%%%% Calculate lossF. 
% load_fTest = chargement(XTest,x_F_min,x_F_max,F,w_max,numTestPoints); %,w_max
% fTest = (1e11*I/(rho*S*w^2))*UxxxxTest-UTest-1/(rho*S*w^2)*load_fTest;% [0,1]  %E*1e11
% 
% zeroTargetTest = zeros(size(fTest), "like", fTest);
% lossTest = mse(fTest, zeroTargetTest);
lossTest = dlarray(0);

% Calculate derivatives U(real) with respect to X
gradientsU = dlgradient(sum(U(1,:),"all"),{X},EnableHigherDerivatives=true);
Ux = gradientsU{1};

gradientsUx = dlgradient(sum(Ux,"all"),{X},EnableHigherDerivatives=true);
Uxx = gradientsUx{1};

gradientsUxx = dlgradient(sum(Uxx,"all"),{X},EnableHigherDerivatives=true);
Uxxx = gradientsUxx{1};

gradientsUxxx = dlgradient(sum(Uxxx,"all"),X,EnableHigherDerivatives=true);
Uxxxx = gradientsUxxx;


% Calculate derivatives U(imaginary) with respect to X
gradientsV = dlgradient(sum(U(2,:),"all"),{X},EnableHigherDerivatives=true);
Vx = gradientsV{1};

gradientsVx = dlgradient(sum(Vx,"all"),{X},EnableHigherDerivatives=true);
Vxx = gradientsVx{1};

gradientsVxx = dlgradient(sum(Vxx,"all"),{X},EnableHigherDerivatives=true);
Vxxx = gradientsVxx{1};

gradientsVxxx = dlgradient(sum(Vxxx,"all"),X,EnableHigherDerivatives=true);
Vxxxx = gradientsVxxx;


% 
% % Calculate derivatives E with respect to X
% 
% gradientsE = dlgradient(sum(E,"all"),{X},EnableHigherDerivatives=true);
% Ex = gradientsE{1};
% 
% gradientsEx = dlgradient(sum(Ex,"all"),X,EnableHigherDerivatives=true);
% Exx = gradientsEx;


% gradientsE = dlgradient(sum(I,"all"),{X},EnableHigherDerivatives=true);
% Ex = gradientsE{1}; 
% gradientsEx = dlgradient(sum(Ex,"all"),X,EnableHigherDerivatives=true);
% Exx = gradientsEx;
% dd_EddU = Exx.*Uxx + 2*Ex.*Uxxx + I.*Uxxxx;
%%%%%% Calculate lossF. 
%load_f = chargement(X,x_F_min,x_F_max,F,w_max,numInternalCollocationPoints); %,w_max

%f =   ((I/(rho*S*w^2))*(Uxxxx.*E*1e11 +2*Ex*1e11.*Uxxx+Exx*1e11.*Uxx))-U-1./(rho*S*w^2).*load_f;% [0,1]  %E*1e11
%f= I*(0.7e11*Uxxxx - 0.7e9*Vxxxx) - U(1,:).*rho*S*w^2;
%I added additional terms such as Vmax Umax to denormalize the data

%%%% getting off gradient error


%real equation
f = I*(Er*1e11*Uxxxx.*Umax - Ei*1e9*Vxxxx.*Vmax) - Umax.*U(1,:).*rho*S*w^2;
zeroTarget = zeros(size(f), "like", f);
lossF = mse(f, zeroTarget);
%imaginary equation
%Ui=U(2,:);
f2=I*(Ei*1e9*Uxxxx.*Umax + Er*1e11*Vxxxx.*Vmax) - Vmax.*U(2,:).*rho*S*w^2;
%f2=I*1e11*(Ei*Uxxxx.*Umax + Er*Vxxxx.*Vmax) - Vmax.*U(2,:).*rho*S*w^2;
zeroTarget2 = zeros(size(f2), "like", f2);
lossF2 = mse(f2, zeroTarget2);

% fbis = (E*Ibis./(rho*Sbis*w^2)).*Uxxxx-U-1./(rho*Sbis*w^2).*load_f;
% zeroTarget = zeros(size(fbis), "like", fbis);
% lossFbis = mse(fbis, zeroTarget);

%lossEregularization = mse(e,ebis);
%%%%% gPINN -> on va dériver le résidu f par rapport aux x_i ,
% effectivement idéalement pour tout x_i, f(x_i) = 0, ainsi df/dx_i = 0
% On ajoute donc cette contrainte.

lossdF = dlarray(0);
if 0
    df = dlgradient(sum(f,"all"),X,EnableHigherDerivatives=true);    
    zeroTarget2 = zeros(size(df), "like", df);
    lossdF = mse(df, zeroTarget2);
end

% Calculate lossObs. Enforce data fitting.

lossObsReal = mse(U_pred_OBS(1,:), ((real((U_OBS)./Umax))));
lossObsImag = mse(U_pred_OBS(2,:), ((imag((U_OBS)./Vmax))));
% lossObsReal = mse(U_pred_OBS(1,:), transpose(stripdims(real((U_OBS)./Umax))));
% lossObsImag = mse(U_pred_OBS(2,:), transpose(stripdims(imag((U_OBS)./Vmax))));
%  lossObsReal = mse(U_pred_OBS(1,:), real(U_OBS));
%  lossObsImag = mse(U_pred_OBS(2,:), imag(U_OBS));

% to test vPINN 

% Combine losses.
%loss = lossObs + lossF + lossU;
loss_tot = 1e-2*lossF ...
          +1e-2*lossF2 ...
        ...+ 1e-3*lossFbis ...
        + 1*lossObsReal ...
        + 1*lossObsImag ...
        + 1e-9*lossdF ...
        ...+ 1e-4*lossEregularization ...
        ...+ 1e-1*lossU_BC0 ...
        ...+ 1e-2*lossdU_BC0 ...
        ...+ 1e-1*lossU_BCL ...
        ...+ 1e-3*lossddU_BCL ...
        ...+ 1e-3*lossdddU_BCL ...
        ;




% Calculate gradients with respect to the learnable parameters.
gradients = dlgradient(loss_tot,parameters);

loss = [(loss_tot)...
    ;(lossF*1e-2)...
    ;(lossF2*1e-2)...
    ...;lossCx... %lossU
    ;(lossObsReal)...
    ;(lossObsImag)...
    ];

end
