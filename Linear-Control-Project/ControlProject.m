M = 0.5;
m = 0.2;
b = 0.1;
I = 0.006;
g = 9.8;
l = 0.3;
q = (M+m)*(I+m*l^2)-(m*l)^2;
A = [0      1              0           0;
     0 -(I+m*l^2)*b/q  (m^2*g*l^2)/q   0;
     0      0              0           1;
     0 -(m*l*b)/q       m*g*l*(M+m)/q  0];
B = [     0;
     (I+m*l^2)/q;
          0;
        m*l/q];
C = [1 0 0 0;
     0 0 1 0];
D = [0;
     0];
s = tf('s');
P_pend = (m*l*s/q)/(s^3 + (b*(I + m*l^2))*s^2/q - ((M + m)*m*g*l)*s/q - b*m*g*l/q);
P_cart = (((I+m*l^2)/q)*s^2 - (m*g*l/q))/(s^4 + (b*(I + m*l^2))*s^3/q - ((M + m)*m*g*l)*s^2/q - b*m*g*l*s/q);
fprintf("poles and zeros of P_pend")
pole(P_pend)
zero(P_pend)
fprintf("poles and zeros of P_cart")
pole(P_cart)
zero(P_cart)
Ap = [0, 1; (m*m*l*l*g)/(I*(M+m)+M*m*l*l), 0];
Bp = [0; (I+m*l*l)/(I*(M+m)+M*m*l*l)];
Cp=[1 0];
Dp=[0];
Poles = [-7.25; -7.5];
k = place(Ap, Bp, Poles)
Acl=Ap-Bp*k;
Ccl=Cp-Dp*k;
syscl= ss(Acl,Bp,Cp,Dp);
stepinfo(syscl)


