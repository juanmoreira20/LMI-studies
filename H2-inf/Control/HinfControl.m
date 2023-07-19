clear all; close all; clc;

A = [2, 1, -2;
     1,-1, -3 ;
     4, 0, -1];
B1 = [1, 0;
    0, 3;
    3, 1];
B2 = [1, 0.2, -0.5]';
C = [2, 1, -0.5];
D1 = [1, 1];
D2 = 0.05;
[aux,n] = size(A);
[aux,m] = size(B1);
[q,aux] = size(C);
gamma = sdpvar(1,1);
X = sdpvar(n,n);
W = sdpvar(m,n);

a11 = (A*X + B1*W)' +A*X + B1*W;
a12 = B2;
a13 = (C*X + D1*W)';
a21 = a12';
a22 = -gamma*eye(q);
a23 = D2';
a31 = a13';
a32 = a23';
a33 = -gamma*eye(q);

LMI = [a11, a12, a13;
       a21, a22, a23;
       a31, a32, a33];
LMI = [LMI <= 0; X>=0];
options = sdpsettings('solver','mosek');
optimize(LMI,gamma,options)
C = checkset(LMI)