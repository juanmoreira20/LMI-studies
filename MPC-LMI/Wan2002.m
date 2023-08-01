close all; clear;clc;
%% Definition of the system in state space
m = 0.140; % Massa da haste
l = 0.135/2; % Tamanho da haste
J = 1/3*m*((0.135/2)^2+(0.018)^2); % Momento de inercia da haste vertical
g = 9.8; % Aceleração da gravidade
R = 4.286; % Resistência do Motor CC
r = 0.156; % Tamanho da haste horizontal
K_m = 0.1797; % Coeficiente de inercia mecânica do motor
K_e = 0.2087; % Coeficiente da FEM do motor
Ts = 0.20; % Tempo de amostragem
I = 1/3*m*((0.135/2)^2+(0.018)^2)+m*(0.135/2)^2; % Momento de inercia da haste horizontal
Q_eq = I*J+I*m*l^2+J*m*r^2; % Quantidade de Energia equivalente
A_22 = -(K_m*K_e*(J+m*l^2))/(Q_eq*R);
A_23 = (m^2*l^2*r*g)/Q_eq;
A_42 = -(m*r*l*K_m*K_e)/(Q_eq*R);
A_43 = (m*g*l*(I+m*r^2))/Q_eq;
B_21 = (K_m*(J+m*l^2))/(Q_eq*R);
B_41 = (m*r*l*K_m)/(Q_eq*R);

% State Space
A = [0 1 0 0; 0 A_22 A_23 0; 0 0 0 1; 0 A_42 A_43 0];
B = [0; B_21; 0; B_41];
D =[0]; C = [1 0 1 0];

% Discretização da planta
sys=ss(A, B, C, D); 
sys=c2d(sys, Ts);
[A,B,C,D] = ssdata(sys);


% Polytopic uncertanty
A1=1.01*A; A2=0.85*A; 

%% Step response
tset=0:1:10; % step response
x0=[0.1745; 0; 0; 0]; % Initial condition
for k=1:1:size(tset,2)-1
    x0(:,k+1)=((A1+A2)/2)*x0(:,k);
end

% Constraints for quadratic supply, control saturation and system decay
uMax = 2;
%Q1 = [50 0 0 0; 0 25 0 0; 0 0 50 0; 0 0 0 10]; not working
Q1 = diag([1 1 1 1]);
R = 1;
p = 0.5; % decay

% LMI formulation
[~,n] = size(A1);
[~,m] = size(B);
Q=sdpvar(n,n,'symmetric');
Q0=sdpvar(n,n,'symmetric');
Y=sdpvar(m,n,'full');
Y0=sdpvar(n,m,'full');
X=sdpvar(m,m,'symmetric');
xk = sdpvar(n,m);
gamma=sdpvar(1,1);

% State Observer:
options = sdpsettings('solver','mosek');
solvesdp([[p*p*Q0 (Q0*A-Y0*C);(Q0*A-Y0*C)' Q0]>=0, Q0>=0],[], options);
Lp=inv(value(Q0))*value(Y0)

% LMI for first polytopic constraint
a11 = Q;
a12 = (A1*Q+B*Y)';
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
% LMI for second polytopic constraint
a12 = (A2*Q+B*Y)';
a21 = a12';

LMI2 = [a11, a12, a13, a14;
        a21, a22, a23, a24;
        a31, a32, a33, a34;
        a41, a42, a42, a44];

LMIs = [LMI1 >=0, LMI2 >=0, Q >= 0, X>=0, Y>=0, gamma >= 0, [1 xk';xk Q]>=0, [X Y;Y' Q]>=0, X<=uMax*uMax];

%% Look up table formulation:
aux1=[]; 
aux2=[]; 
Table=[];
for k=1:1:size(x0,2)
    model = optimizer(LMIs, gamma,options,xk,{X,Y,Q});
    x=[x0(1,k) x0(2,k) x0(3,k) x0(4,k)]';
    QY = model{x}; % solver for each line
    
    QQ(:,:,k)=QY{3}; YY(:,:,k)=QY{2};
    F(:,:,k) = YY(:,:,k)*inv(QQ(:,:,k));
    
    %% Lookup-table
    Table=[Table; F(:,:,k)];
    
    aux1=[aux1; F(1,1,k) ];
    aux2=[aux2; F(1,2,k) ];
    LMIs=[LMIs, Q<=QQ(:,:,k)]; % acrescimo da restrição do algoritmo offline
end

%% Robustness
Qpoly=sdpvar(2*n, 2*n,'symmetric');

Poly = Qpoly>=0; 
for i=1:1:length(F)
    %% Augmented system (with observed states)
    Apoly1=[A1 B*F(:,:,k); Lp*C A1+B*F(:,:,i)-Lp*C];
    Apoly2=[A2 B*F(:,:,k); Lp*C A2+B*F(:,:,i)-Lp*C];
    
    %% LMI formulation
    eqpoly1 = [Qpoly (Qpoly*Apoly1)';(Qpoly*Apoly1) Qpoly];
    eqpoly2 = [Qpoly (Qpoly*Apoly2)';(Qpoly*Apoly2) Qpoly];
    
    Poly =[Poly, eqpoly1>=0, eqpoly2>=0];
end
solvesdp(Poly,[],options), % Uso do solver

%% Time formulation
ts=0.1;  
t=0:ts:30; 
% Look-up table line
N= randi([1,10]);      
% Initial conditions
x=[0; 0; 0.1745; 0]; 
% x estimated
xest=[0; 0; 0; 0];
for k=1:1:length(t)-1
    %% V1
    u(k)=F(:,:,N)*xest(:,k);
    xest(:,k+1)=A*xest(:,k)+B*u(k)+Lp*C*(x(:,k)-xest(:,k));
    x(:,k+1)=A*x(:,k)+B*u(k);
    
end

%% Estados
figure

plot(t, x(1,:),'r','linewidth',2);
hold on
plot(t, x(2,:),'b','linewidth',2); 
hold on
plot(t, x(3,:),'g','linewidth',2); 
hold on
plot(t, x(4,:),'k','linewidth',2); 
hold on
grid on
