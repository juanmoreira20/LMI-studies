%Example of D-stability using a area of D(q,r) ∩ H(alpha) such as
% D(q,r) = {x+yj / (x+q)^2 + y^2 < r^2}
% H(alpha) = {x+yj / x < -alpha}
% ẋ = Ax(t)
% with those regions we must use diffent constrains, such as
% for H: (A + αI )'P + P(A + αI) < 0
% for D: (A + qI )'P(A + qI ) − r^2P < 0
A = [-4.2386 -0.2026 0.7193
2.6649 -2.8342 0.0175
0.0344 0.0005 -3.1772];
[~,n]=size(A);
P = sdpvar(n,n);
options = sdpsettings('solver','sedumi');
alpha = 2;
q = 3;
r = 3;
LMI = [P>=0,((A + alpha*eye(n,n) )'*P) + (P*(A + alpha*eye(n,n))) <=0 ,(A + q*eye(n,n) )'*P*(A + q*eye(n,n) ) - (r^2)*P <= 0];
optimize(LMI,[],options);
value(P)
eig(A)