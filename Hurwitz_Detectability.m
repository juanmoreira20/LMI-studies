%Hurwitz Detectability: A linear(A,B,C) system is hurwitz stabilizable if
%and only if its dual counterpart (A',B',C') is Hurwitz Detectable
A = [-0.0558 -0.9968 0.0802 0.0415
0.5980 -0.1150 -0.0318 0.0000
-3.0500 0.3880 -0.4650 0.0000
0.0000 0.0805 1.0000 0.0000];
C = [0 1 0 0
0 0 0 1];
options = sdpsettings('solver','sedumi');
[~,n]=size(A);
P = sdpvar(n,n);
LMI = [(A'*P)+(P*A) <= C'*C, P>=0];
optimize(LMI,[],options);
value(P)