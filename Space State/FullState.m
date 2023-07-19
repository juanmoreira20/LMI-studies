clear;
close all;
clc;

% Xb. = A*Xb + B*u + L(yb - y)
% Xb. = A*Xb + B*u + L*C*Xb - L*C*X
% X.  = A*X +  B*u 
% e = X - Xb => e. = (A+L*C)*e, logo, para que o erro tenda a 0, os
% autovalores de A+L*C devem estar no semiplano esquerdo, logo, podemos
% fazer uma analogia com A+B*K, tal que as LMI's resolvem o problema:
% LMI = [P*A+A'*P+W*'+C'*W'<=0, P>=0]
% L = inv(P)*W;
A = [-0.5 0 0;
      0 -2 10;
      0  1 -2];
B = [1 0; -2 2; 0 1];
C = [1 0 0; 0 0 1];
[~,n] = size(A);
[~,m]= size(B);
P = sdpvar(n,n,'symmetric');
W = sdpvar(n,m);
LMI = [P*A+A'*P+W*C+C'*W'<=0, P>=0];
options = sdpsettings('solver','mosek');
optimize(LMI,[],options)
L = inv(value(P))*value(W);
eig(A+L*C)
    