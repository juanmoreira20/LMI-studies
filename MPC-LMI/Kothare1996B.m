close all; clear;clc;
%% Definition of the system in state space
% physical formulation
k = 0.787;
% polytopic uncertanties
alphaMin = 0.1;
alphaMax = 10;
A1 = [1 0.1;
      0 1-0.1*alphaMin];
A2 = [1 0.1;
      0 1-0.1*alphaMax];

% state space formulation
B = [0, 0.1*k]';
C = [1 0];
D = 0;

% Constraints for quadratic supply rate and control saturation
uMax = 2;
Q1 = [1,0;
      0,0];
R = 0.00002;

% LMI formulation
[~,n] = size(A1);
[~,m] = size(B);
Q=sdpvar(n,n,'symmetric');
Y=sdpvar(m,n,'full');
X=sdpvar(m,m,'symmetric');
xk = sdpvar(n,m);
gamma=sdpvar(1,1);

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
LMI = [LMI1 >=0, LMI2 >=0, Q >= 0, X>=0, gamma >= 0, [1 xk';xk Q]>=0, [X Y;Y' Q]>=0, X<=uMax*uMax] 
% Optimizer formulation
options = sdpsettings('solver','mosek');
model = optimizer(LMI, gamma,options,xk,{Y,Q});

%% Initial Conditions
x = [1; 0];
xs = x;
QY = model{x};
F = QY{1}*inv(QY{2});
Fs = F;

tsim = 10; ts=0.1;
t = 0:ts:(tsim);

F'

for k=1:ceil(tsim/ts)
    alpha(k)=10*rand; % alfa randomico
    A = [1 0.1; 0 1-0.1*alpha(k)];
    u(k) =  F*x(:,k);
    x(:,k+1) = A*x(:,k) + B*u(k);
% Online formulation:
    QY = model{x(:,k+1)};
    F = QY{1}*inv(QY{2});
    
    us(k) = Fs*(xs(:,k));
    xs(:,k+1) = A*xs(:,k) + B*us(k);
end

%%
%% Plot
t = 0:ts:(tsim);

figure();
subplot(1,2,1);
plot(t, x(1,:),'k', t,xs(1,:),'k--');
xlabel('Time (sec)'); ylabel('\theta (rad)');
legend('MPC Online','MPC Offline');

subplot(1,2,2);
plot(t(1:end-1),u,'k-',t(1:end-1),us,'k--');
xlabel('Time (sec)'); ylabel('u (volts)');
legend('MPC Online','MPC Offline');
