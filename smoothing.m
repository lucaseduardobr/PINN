function [x,y] = smoothing(X,U,n)
i=0;
b=0;
[X,U]=fixOri(X,U);

while b<=0.99
    i=i+1;
    f_n="fourier"+num2str(i);
    try
    [f,gof] = fit(X',U',f_n);
    catch
    [f,gof] = fit(X,U,f_n);    
    end
    b=gof.adjrsquare;

end

x_min=min(X);
x_max=max(X);
x=linspace(x_min,x_max,n);
y=feval(f,x)';
   

end