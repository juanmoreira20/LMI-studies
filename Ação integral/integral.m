clear;
clc;
close all;
H = tf(1,[1 0.1 2]);
[num,den] = tfdata(H,'v');
[A,B,C,D] = tf2ss(num,den);
Aa = [A,zeros(size(A,1),size(C,1));-C,zeros(size(C,1),size(C,1))];
Bu = [B;zeros(size(A,1)+size(C,1)-size(B,1),size(B,2))];
Ca = [C 0];

[~,n] = size(Aa);
[~, m] = size(Bu);

Q = sdpvar(n,n,'symmetric');
N = sdpvar(m,n);

LMI =[Q>=0;
    Q*Aa'+N'*Bu'+Aa*Q+Bu*N <=0];
%resolvendo LMI
options = sdpsettings('solver','sedumi');
optimize(LMI,[],options)
Qv = value(Q);
Nv = value(N);
%ganho aumentado
Ka =  Nv*inv(Qv);

% ganho malha fechada
K = Ka(1:size(A,1));
% ganho integral
H = Ka(size(A,1)+1:end);
Tsim = 40;
