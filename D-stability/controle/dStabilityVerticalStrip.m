clear;
clc;
close all;
H = tf(1,[1 0.1 2]);
[num,den] = tfdata(H,'v');
[A,B,C,D] = tf2ss(num,den);
a = -1;
b = -2;
[~,n] = size(A);
[~, m] = size(B);

P = sdpvar(n,n,'symmetric');
W = sdpvar(m,n);
L1 = [A'*P + P*A +B*W + W'*B' - 2*a*P];
L2 = -[A'*P + P*A +B*W + W'*B' - 2*b*P];

LMI = [P>=0;
        L1 <=0;
        L2 <=0];
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

