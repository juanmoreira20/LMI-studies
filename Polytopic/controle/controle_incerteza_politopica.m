clear;
clc;
close all;
%%descrição do sistema
%matrizes de teste
At0 = [-4 4;-5 0];
At1 = [-4 4;-5 0] + [-2 2;-1 4];
At2 = [-4 4;-5 0] - [-2 2;-1 4];
B = [1;2];
C = [1 1];
D = 0;
%incerteza pertence a (-1,1)
inc = 1 - 2*rand; 
% matriz da incerteza principal
A = [-4 4;-5 0] + inc*[-2 2;-1 4];
%LMI infeasible, logo, devemos utilizar alguma estratégia de controle para
%estabilizar dentre a região indicada


%% Estratégia de controle: Se podemos testar a estabilidade por vértice, também podemos criar a lei de controle através deles
%Utilizando a estratégia de controle contínuo por LMI temos:
%Espaço no Rn
[~,n] = size(A);
[~, m] = size(B);

% variáveis de decisão
P = sdpvar(n,n);
W = sdpvar(m,n);

%Sintese do controle
LMI = [ P>=0 ;
    P*At1'+At1*P+W'*B'+B*W<=0;
    P*At2'+At2*P+W'*B'+B*W<=0];
options = sdpsettings('solver','sedumi');
optimize(LMI,[],options)
Pv = value(P);
Wv = value(W);
K = Wv*inv(Pv);
checkset(LMI)
% sinal de saida e controle incerteza
Outinc = ss((A+B*K),B,C,D); 
Csinc = ss((A+B*K),B,K,D); 

% sinal de saida e controle incerteza
Out1 = ss((At1+B*K),B,C,D); 
Cs1 = ss((At1+B*K),B,K,D); 

% sinal de saida e controle incerteza
Out2 = ss((At2+B*K),B,C,D); 
Cs2 = ss((At2+B*K),B,K,D); 

figure
title('Malha fechada A')
step(Outinc,'r')
hold on
step(Csinc,'y')
legend('Saída','Controle')

figure
title('Malha fechada At1')
step(Out1,'r')
hold on
step(Cs1,'y')
legend('Saída','Controle')

figure
title('Malha fechada At2')
step(Out2,'r')
hold on
step(Cs2,'y')
legend('Saída','Controle')