clear;
close all;
clc;

%processo
Ts = 0.1;

H = tf(1,[1 -1 2]);
Hd = c2d(H,Ts);

[num,den] = tfdata(Hd,'v');

[A,B,C,D] = tf2ss(num,den);

%Espaço no Rn
[~,n] = size(A);
[~, m] = size(B);

% variáveis de decisão
P = sdpvar(n,n);
W = sdpvar(m,n);

%sintese do controle pós aplicação do complemento de schur 
LMI = [ P>=0 ; 
    [-P, P*A'+W'*B';A*P+B*W, -P] <= 0];
%%lei de controle simplificada
% LMI = [ P>=0 ; 
%     [-P, P*A';A*P+B*W, -P-B*B'] <= 0];

%resolvendo a lmi
options = sdpsettings('solver','sedumi');
optimize(LMI,[],options)

%Setando ganho
Pv = value(P);
Wv = value(W);
K = Wv/Pv;
%%ganho da lei simplificada
%K = −inv(2*eye(n,n) + B'*inv(P)*B)*B'*inv(P)*A

%sinal de saída e controle em espaço de estados
Out = ss((A+B*K),B,C,D,Ts); 
CS = ss((A+B*K),B,K,0,Ts); 

%figuras
figure
step(Out,'r')
hold on
step(CS,'y')
legend('saida','sinal de controle')

figure
step(Hd)

