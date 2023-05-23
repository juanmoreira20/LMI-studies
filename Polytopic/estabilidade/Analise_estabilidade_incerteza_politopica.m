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
% Análise de estabilidade do sistema
%% Para o sistema sem incerteza:
%Espaço no Rn
[~,n] = size(A);
[~, m] = size(B);

% variáveis de decisão
P = sdpvar(n,n);
LMI = [P >=0;
        At0*P+P*At0' <=0];
options = sdpsettings('solver','sedumi');
optimize(LMI,[],options)
Pv = value(P);
checkset(LMI);
eig(At0);

%% Para o sistema com incerteza igual a 1
% variáveis de decisão
P = sdpvar(n,n);
LMI = [P >=0;
        At1*P+P*At1' <=0];
options = sdpsettings('solver','sedumi');
optimize(LMI,[],options)
Pv = value(P);
checkset(LMI);
eig(At1);

%% Para o sistema com incerteza igual a -1
% variáveis de decisão
P = sdpvar(n,n);
LMI = [P >=0;
        At2*P+P*At2' <=0];
options = sdpsettings('solver','sedumi');
optimize(LMI,[],options)
Pv = value(P);
checkset(LMI);
eig(At2);
%% Para incerteza de -1 a 1(Apenas necessário testar os vértices):
P = sdpvar(n,n);
LMI = [P >=0;
        At2*P+P*At2' <=0;
        At1*P+P*At1' <=0];
options = sdpsettings('solver','sedumi');
optimize(LMI,[],options)
Pv = value(P);
checkset(LMI);
%LMI infeasible, logo, devemos utilizar alguma estratégia de controle para
%estabilizar dentre a região indicada

% sinal de saida incerteza
Outinc = ss(A,B,C,D); 
 

% sinal de saida  incerteza
Out1 = ss(At1,B,C,D); 


% sinal de saida  incerteza
Out2 = ss(At2,B,C,D); 


figure
title('Malha aberta A')
step(Outinc,'r')

figure
title('Malha aberta At1')
step(Out1,'r')


figure
title('Malha aberta At2')
step(Out2,'r')






