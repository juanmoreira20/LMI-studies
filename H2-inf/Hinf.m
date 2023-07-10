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
gamma = sdpvar(1,1);
P = sdpvar(n,n);

a11 = A'*P+P*A;
a12 = P*B;
a13 = C';
a21 = B'*P;
a22 = -gamma*eye(q);
a23 = D';
a31 = C;
a32 = D;
a33 = -gamma*eye(q);

LMI = [a11, a12, a13;
       a21, a22, a23;
       a31, a32, a33];
LMI = [LMI <= 0; P>=0];
options = sdpsettings('solver','mosek');
optimize(LMI,gamma,options)
C = checkset(LMI)