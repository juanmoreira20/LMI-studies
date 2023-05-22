% Hurwitz Stabilizability to find a control law(gain) for the proposed
% closed loop of a linear system:
% ẋ = Ax(t)+Bu(t)
% y = Cx(t)
A = [-1 0 1 1 0 0 -1 2
0 0 1 0 -3 -1 0 0
-1 -2 1 0 -3 -1 0 0
0 0 0 0 1 0 0 0
0 0 -1 0 -3 0 4 2
0 0 0 0 0 0 1 0
0 0 0 0 0 0 0 1
-1 2 0 0 1 -1 -1 2];
B = [1 0 0 0
0 0 0 0
0 1 0 0
0 0 0 0
0 0 1 0
0 0 0 0
0 0 0 0
0 0 0 0];
options = sdpsettings('solver','sedumi');
[~,n]=size(A);
P = sdpvar(n,n);
LMI = [(A*P)+(P*A') <= B*B', P>=0];
optimize(LMI,[],options);
value(P)
K = -(1/2)*B'*inv(value(P));
eig(A+B*K)

