clc;
clear all;
close all;
syms wn;
syms zta;
wn = 0.8*(zta)^-1;
theta = acos(zta);
b = (1 - (zta)^2)^0.5;
eq1 = (pi - theta)/(wn*b) - 0.5;
ztasol = solve(eq1,zta);
wn = 0.8*(ztasol)^-1;
syms s;
eq2 = s^2 + 2*ztasol*wn*s + (wn)^2;
ssolv = solve(eq2,s);
%% state feedback
%state matrix
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
C = [1 0 0 0;
     0 0 1 0];
D = [0;
     0];
% feedback methode
pm = [-0.8-3.5811i -0.8+3.5811i -10 -8];
k = place(A,B,pm);
A = A - B*k ;


[anew , bnew] = ss2tf(A ,B ,C ,D);
