clc;
clear all;
close all;
x = tf([4.182e-06 0 - 0.0001025],[2.3e-06 4.182e-07 -7.172e-05 -1.025e-05 0])
phi = tf([4.545 -1.277e-16],[1 0.1818 -31.18 -4.455])
px = pole(x)
pphi = pole(phi)
zx = zero(x)
zp = zero(phi)
subplot(2 , 2 ,1);
pzmap(x);
subplot(2 , 2 ,2);
pzmap(phi);
subplot(2 , 2 ,3);
bode(x);
subplot(2 , 2 ,4);
bode(phi);
%% define the conditions
syms x
eqn = ((pi - (1 - x)^0.5)/((0.8/x)*(acos(x))) - 0.5);

% Solve the equation
sol = vpasolve(eqn, x);

% Display the solutions
disp(sol);
%% state feedback
M = .5;
m = 0.2;
b = 0.1;
I = 0.006;
g = 9.8;
l = 0.3;

p = I*(M+m)+M*m*l^2; %denominator for the A and B matrices

A = [0      1              0           0;
     0 -(I+m*l^2)*b/p  (m^2*g*l^2)/p   0;
     0      0              0           1;
     0 -(m*l*b)/p       m*g*l*(M+m)/p  0];
B = [     0;
     (I+m*l^2)/p;
          0;
        m*l/p];
ptwo = [- 0.8 - 3.3052267077375095057033673755362i ,- 0.8 + 3.3052267077375095057033673755362i, -8 ,-10 ];
K = place(A ,B ,ptwo);
A = A - B*K;
C = [1 0 0 0;
     0 0 1 0];
D = [0;
     0];

[anew , bnew] = ss2tf(A ,B ,C ,D);
