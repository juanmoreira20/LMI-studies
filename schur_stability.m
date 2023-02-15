%stability of a discrete x(k+1) = Ax(k) system using Lyapunov approach of
%schur criteria
A = [-0.5 2 0
    0 -0.25 0
    1 0 -0.5];
options = sdpsettings('solver','sedumi');

[~,n]=size(A);
P = sdpvar(n,n);
LMI = [A*P*A' - P <= 0, P>=0];
optimize(LMI,[],options);
value(P)
eig(A)
