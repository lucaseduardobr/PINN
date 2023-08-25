syms x a0 a1 a2  b1 b2  w

v_k = -(a0 + a1*cos(x.*w) + b1*sin(x.*w) + a2*cos(2*x.*w) + b2*sin(2*x.*w));

v1=diff(v_k,x,1);

v1=a1*w*sin(w*x) - 2*b2*w*cos(2*w*x) - b1*w*cos(w*x) + 2*a2*w*sin(2*w*x);


v2=diff(v_k,x,2);

v2=a1*w^2*cos(w*x) + 4*a2*w^2*cos(2*w*x) + b1*w^2*sin(w*x) + 4*b2*w^2*sin(2*w*x);

v3=diff(v_k,x,3);

v3=b1*w^3*cos(w*x) + 8*b2*w^3*cos(2*w*x) - a1*w^3*sin(w*x) - 8*a2*w^3*sin(2*w*x);

v4=diff(v_k,x,4);

v4=- a1*w^4*cos(w*x) - 16*a2*w^4*cos(2*w*x) - b1*w^4*sin(w*x) - 16*b2*w^4*sin(2*w*x);