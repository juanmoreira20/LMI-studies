clear;
close all;
clc;
%system
A = [-0.5 0 0;
      0 -2 10;
      0  1 -2];
B = [1 0; -2 2; 0 1];
C = [1 0 0; 0 0 1];
A12 = [0;1];
A22 = -2;
[n,~] = size(A12);
[~,m]= size(A22);
P = sdpvar(m,m,'symmetric');
W = sdpvar(m,n);
LMI = [P*A22+A22'*P+W*A12+A12'*W'<=0, P>=0];
options = sdpsettings('solver','mosek');
optimize(LMI,[],options)
L = inv(value(P))*value(W);
eig(A+L*C)