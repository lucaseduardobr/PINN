function Uxxxx = deriveeOrdre4 (U,X)
nx = length(X);
N = nx;
dx = 1/(nx-1);

%num = [-2,-1, 0, 1, 2];
num = [-3, -2, -1, 0, 1, 2, 3];
A=spalloc(N,N,7*N);
%A = zeros(N,N);
for k=4:nx-3
        dx = X(k+1)-X(k);
        %coeff=[1/dx^4,-4/dx^4,6/dx^4,-4/dx^4,1/dx^4];
        coeff = 1/(6*dx^4)*[-1 12 -39 56 -39 12 -1];
        A(k,k+num)=coeff;
end
A;
Uxxxx = A*U';
end