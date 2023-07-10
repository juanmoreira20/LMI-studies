clear all; close all; clc;
% qsi = 0.1;
% w = 1;
% k = 1;
qsi = 1;
w = 3;
k = 20;

A = [-2*qsi*w, -w*w;
     1,          0];
B = [1, 0]';
C = [0, k];
[aux,n] = size(A);
[aux,m] = size(B);
[q,aux] = size(C);
D = 0*sdpvar(q,m);
rho = sdpvar(1,1);
X = sdpvar(n,n);

LMI = A*X+X*A'+B*B';

LMI = [LMI <= 0; X>=0, trace(C*X*C') <= rho];
options = sdpsettings('solver','mosek');
optimize(LMI,rho,options)
C = checkset(LMI)