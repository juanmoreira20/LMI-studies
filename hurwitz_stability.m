%stability of a continuous ẋ = Ax(t) system using Lyapunov approach of
%Huritz criteria
A = [-0.5 1 0 0 0
    -62 -0.16 0 0 30
    10 0 -50 40 5
    0 0 0 -40 0
    0 0 0 0 -20];
options = sdpsettings('solver','sedumi');
[~,n]=size(A);
P = sdpvar(n,n);
LMI = [(A'*P)+(P*A) <= 0, P>=0];
optimize(LMI,[],options);
value(P)
eig(A)