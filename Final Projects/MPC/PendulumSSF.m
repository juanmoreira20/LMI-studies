close all; clear;clc;
tic
%% Definition of the system in state space
% Simulação do modelo do pendulo invertido Rotativo
m = 0.140; % Massa da haste
m = m -0.7*m;
l = 0.135/2; % Tamanho da haste
J = 1/3*m*((0.135/2)^2+(0.018)^2); % Momento de inercia da haste vertical
g = 9.8; % aceleração da gravidade
R = 4.286; % Resistência do Motor CC
r = 0.156; % tamanho da haste horizontal
K_m = 0.1797; % Coeficiente de inercia mecânica do motor
K_e = 0.2087; % Coeficiente da FEM do motor
I = 1/3*m*((0.135/2)^2+(0.018)^2)+m*(0.135/2)^2;% Momento de inercia da haste horizontal
Q_eq = I*J+I*m*l^2+J*m*r^2; %Quantidade de Energia equivalente
A_22 = -(K_m*K_e*(J+m*l^2))/(Q_eq*R);
A_23 = (m^2*l^2*r*g)/Q_eq;
A_42 = -(m*r*l*K_m*K_e)/(Q_eq*R);
A_43 = (m*g*l*(I+m*r^2))/Q_eq;
B_21 = (K_m*(J+m*l^2))/(Q_eq*R);
B_41 = (m*r*l*K_m)/(Q_eq*R);
% Modelo no espaço de estado
A = [0 1 0 0; 0 A_22 A_23 0; 0 0 0 1; 0 A_42 A_43 0];
B = [0; B_21; 0; B_41];
C = [1 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 0];
D =[0 0 0 0]';
ts = 0.02;
sys = ss(A,B,C,D);
[A,B,C,D] = ssdata(c2d(sys,ts));

% Constraints for quadratic supply rate and control saturation
uMax = 3;
Q1 = [50 0 0 0; 0 25 0 0; 0 0 50 0; 0 0 0 10]; %Matriz de Ponderação
%Q1 = eye(4);
R = 1;
% LMI formulation
[~,n] = size(A);
[~,m] = size(B);
Q=sdpvar(n,n,'symmetric');
Y=sdpvar(m,n,'full');
X=sdpvar(m,m,'symmetric');
xk = sdpvar(n,m);
gamma=sdpvar(1,1);

% LMI for first polytopic constraint
a11 = Q;
a12 = (A*Q+B*Y)';
a13 = (sqrt(Q1)*Q)';
a14 = (sqrt(R)*Y)';
a21 = a12';
a22 = Q;
a23 = 0*Q;
a24 = 0*Y';
a31 = a13';
a32 = a23';
a33 = gamma*eye(n);
a34 = 0*Y';
a41 = a14';
a42 = a24';
a43 = a34';
a44 = gamma*eye(m);

LMI1 = [a11, a12, a13, a14;
        a21, a22, a23, a24;
        a31, a32, a33, a34;
        a41, a42, a42, a44];
m = 0.14;
m = m +0.7*m;
l = 0.135/2; % Tamanho da haste
J = 1/3*m*((0.135/2)^2+(0.018)^2); % Momento de inercia da haste vertical
g = 9.8; % aceleração da gravidade
R = 4.286; % Resistência do Motor CC
r = 0.156; % tamanho da haste horizontal
K_m = 0.1797; % Coeficiente de inercia mecânica do motor
K_e = 0.2087; % Coeficiente da FEM do motor
I = 1/3*m*((0.135/2)^2+(0.018)^2)+m*(0.135/2)^2;% Momento de inercia da haste horizontal
Q_eq = I*J+I*m*l^2+J*m*r^2; %Quantidade de Energia equivalente
A_22 = -(K_m*K_e*(J+m*l^2))/(Q_eq*R);
A_23 = (m^2*l^2*r*g)/Q_eq;
A_42 = -(m*r*l*K_m*K_e)/(Q_eq*R);
A_43 = (m*g*l*(I+m*r^2))/Q_eq;
B_21 = (K_m*(J+m*l^2))/(Q_eq*R);
B_41 = (m*r*l*K_m)/(Q_eq*R);

A = [0 1 0 0; 0 A_22 A_23 0; 0 0 0 1; 0 A_42 A_43 0];
B = [0; B_21; 0; B_41];
C = [1 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 0];
D =[0 0 0 0]';
ts = 0.02;
sys = ss(A,B,C,D);
[A,B,C,D] = ssdata(c2d(sys,ts));
[~,m] = size(B);
a11 = Q;
a12 = (A*Q+B*Y)';
a13 = (sqrt(Q1)*Q)';
a14 = (sqrt(R)*Y)';
a21 = a12';
a22 = Q;
a23 = 0*Q;
a24 = 0*Y';
a31 = a13';
a32 = a23';
a33 = gamma*eye(n);
a34 = 0*Y';
a41 = a14';
a42 = a24';
a43 = a34';
a44 = gamma*eye(m);


LMI2 = [a11, a12, a13, a14;
        a21, a22, a23, a24;
        a31, a32, a33, a34;
        a41, a42, a42, a44];
LMI = [LMI1 >=0,LMI2>=0, Q >= 0, X>=0, gamma >= 0, [1 xk';xk Q]>=0, [X Y;Y' Q]>=0, X <= uMax^2] ;
% Optimizer formulation
options = sdpsettings('solver','mosek');
model = optimizer(LMI, gamma,options,xk,{Y,Q});

%% Condições Iniciais
x = [0,0,0.1745,0]';
QY = model{x};
F = QY{1}*inv(QY{2});

tsim = 7;
t = 0:ts:(tsim);

F'
m = 0.14; % Massa da haste
l = 0.135/2; % Tamanho da haste
J = 1/3*m*((0.135/2)^2+(0.018)^2); % Momento de inercia da haste vertical
g = 9.8; % aceleração da gravidade
R = 4.286; % Resistência do Motor CC
r = 0.156; % tamanho da haste horizontal
K_m = 0.1797; % Coeficiente de inercia mecânica do motor
K_e = 0.2087; % Coeficiente da FEM do motor
I = 1/3*m*((0.135/2)^2+(0.018)^2)+m*(0.135/2)^2;% Momento de inercia da haste horizontal
Q_eq = I*J+I*m*l^2+J*m*r^2; %Quantidade de Energia equivalente
A_22 = -(K_m*K_e*(J+m*l^2))/(Q_eq*R);
A_23 = (m^2*l^2*r*g)/Q_eq;
A_42 = -(m*r*l*K_m*K_e)/(Q_eq*R);
A_43 = (m*g*l*(I+m*r^2))/Q_eq;
B_21 = (K_m*(J+m*l^2))/(Q_eq*R);
B_41 = (m*r*l*K_m)/(Q_eq*R);
% Modelo no espaço de estado
A = [0 1 0 0; 0 A_22 A_23 0; 0 0 0 1; 0 A_42 A_43 0];
B = [0; B_21; 0; B_41];
C = [1 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 0];
D =[0 0 0 0]';
ts = 0.02;
sys = ss(A,B,C,D);
[A,B,C,D] = ssdata(c2d(sys,ts));
for k=1:ceil(tsim/ts)
    u(k) =  F*x(:,k);
    x(:,k+1) = A*x(:,k) + B*u(k);
% Online formulation:
    QY = model{x(:,k+1)};
    F = QY{1}*inv(QY{2});
end

%%
figure();
plot(t, x(1,:), 'b', t, x(2,:), 'r',t, x(3,:), 'k', t, x(4,:), '--k');
%text(0.5,0.05,'\theta'); text(0.5,-0.25,'d\theta/dt');
ylabel('\theta (rad) and  d\theta/dt (rad/sec)');
xlabel('time(sec)');

toc