clear;
close all;
clc;

%processo

H = tf(1,[1 0.1 2]); % process
[num,den] = tfdata(H,'v');

[A,B,C,D] = tf2ss(num,den);


%Espaço no Rn
[~,n] = size(A);
[~, m] = size(B);

% variáveis de decisão
P = sdpvar(n,n);
W = sdpvar(m,n);

%Sintese do controle
LMI = [ P>=0 ;
    P*A'+A*P+W'*B'+B*W<=0];

%% Lei de controle simplificada
% LMI = [ P>=0 ;
%     P*A'+A*P - B*B'<=0];

%resolvendo LMI
options = sdpsettings('solver','sedumi');
optimize(LMI,[],options)
Pv = value(P);
Wv = value(W);
%ganho
K =  Wv*inv(Pv);
%% ganho simplificado
% K = -(1/2)*B'*inv(Pv);
% sinal de saida e controle
Out = ss((A+B*K),B,C,D); 
Cs = ss((A+B*K),B,K,D); 

figure
step(Out,'r')
hold on
step(Cs,'y')
legend('Saída','Controle')

figure
step(H)