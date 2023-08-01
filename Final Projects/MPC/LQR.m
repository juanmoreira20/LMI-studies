clc, clear all, clc
% Simulação do modelo do pendulo invertido Rotativo
m = 0.140; % Massa da haste
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
Ts = 0.020; % Tempo de amostragem
t = 0:Ts:7; %intervalo de tempo
u = zeros(size(t)); %valor do sinal de controle
[G,H] = c2d(A,B,Ts); % Discretização ZOH Equivalente
x0 = [0; 0; 0.1745; 0]; % Condição inicial a alpha=10º
Tc = ctrb(G,H); %Controlabilidade
if (rank(Tc)==4)
fprintf('O sistema é controlável? \n ');
Q = [50 0 0 0; 0 25 0 0; 0 0 50 0; 0 0 0 10]; %Matriz de Ponderação
R = 1;
K = dlqr(G,H,Q,R); % Ganho LQR discreto
G2 = G-H*K;
y = dlsim(G2,H,C,D,u,x0);
figure(1)
% subplot(3,1,1)
% plot(t,y(:,1),'b-','linewidth',2);
% ylabel('\theta(rad)');
% grid on
% subplot(3,1,2)
% plot(t,y(:,3),'b-','linewidth',2);
% ylabel('\alpha(rad)');
% grid on
% subplot(3,1,3)
% plot(t,-K*y','b-','linewidth',2);
% xlabel('t(s)');
% ylabel('u_a (V)');
% grid on
plot(t, y(1,:), 'b', t, y(2,:), 'r',t, y(3,:), 'k', t, y(4,:), '--k');
%text(0.5,0.05,'\theta'); text(0.5,-0.25,'d\theta/dt');
ylabel('\theta (rad) and  d\theta/dt (rad/sec)');
xlabel('time(sec)');
end