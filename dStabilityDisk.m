clear;
clc;
close all;
H = tf(1,[1 0.1 2]);

r = 0.5;
q = 1;

[num,den] = tfdata(H,'v');
[A,B,C,D] = tf2ss(num,den);

[~,n] = size(A);
[~, m] = size(B);

P = sdpvar(n,n,'symmetric');
W = sdpvar(m,n);
A11 = -r*P;
A21 = q*P+P*A'+W'*B';
A12 = A21';
A22 = A11;
M = [A11 A12; A21 A22];
LMI = [P>=0;
    M<=0];
options = sdpsettings('solver','sedumi');
optimize(LMI,[],options)
Pv = value(P);
Wv = value(W);

K = Wv*inv(Pv);

Out = ss((A+B*K),B,C,D); 
Cs = ss((A+B*K),B,K,D); 

figure
step(Out,'r')
hold on
step(Cs,'y')
legend('Saída','Controle')

figure
step(H)

