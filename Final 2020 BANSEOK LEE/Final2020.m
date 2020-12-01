% Fianl BANSEOK LEE 605351891
clear all
close all

% Define Constants
n = 5; 
m = 0.001; % kg 
l = 0.1; % m
cni = 0.009; % N/(m/s)
cti = 0.0006; % N/(m/s)
ct0 = 0.0006; % N/(m/s)
li = l/(2*n+2); mi = m/(n+1); 
l0 = l/(2*n+2); m0 = m/(n+1);
v0 = -0.15; % m/s

% Define matrices
B = diag(-ones(n,1))+diag(ones(n-1,1),1);
A = diag(ones(n,1))+diag(ones(n-1,1),-1);
e = ones(n,1); 
L = diag(li*ones(n,1)); M = diag(mi*ones(n,1));
Ct = diag(cti*ones(n,1)); Cn = diag(cni*ones(n,1));
C0 = Cn-Ct; 
F = inv(B')*A*L;
h = F'*M*e; 
Lamda = F'*C0 + diag((F'-e*h'/m)*Ct*e);
J = L*M*L/3 + F'*M*F;
D = L*Cn*L/3 + F'*Cn*F;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AA = [zeros(n) eye(n);-B'*inv(J)*v0*Lamda*inv(B') -B'*inv(J)*D*inv(B')];
BB = [zeros(n) ; B'*inv(J)*B];
CC = [eye(n) zeros(n)];
P = [1 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 1 0 0 0 0;
    0 1 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 1 0 0 0;
    0 0 1 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 1 0 0;
    0 0 0 1 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 1 0;
    0 0 0 0 1 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 1;];
AAA = P*AA*inv(P);
BBB = P*BB;
CCC = CC*inv(P);
load('finaldata.mat')
CheckAp = norm(AAA-Ap)/norm(Ap);
CheckBp = norm(BBB-Bp)/norm(Bp);
CheckCp = norm(CCC-Cp)/norm(Cp);
CheckAp
CheckBp
CheckCp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gamma = 0.5*[-1i 1; 1i 1];
xhat = P*[eye(5);1i*w*eye(5)]*phih;
uhat = inv(B)*(-J*w^2+1i*w*D+v0*Lamda)*inv(B')*phih;
X = [xhat conj(xhat)]*Gamma;
U = [uhat conj(uhat)]*Gamma;
O = inv(Gamma)*[1i*w 0;0 -1i*w]*Gamma;
Checkreg = norm(Ap*X+Bp*U-X*O);
Checkreg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options=odeset('RelTol',1e-4,'Refine',5);
simtime=[0 2];
xinit=[1;1;0;0;0;0;0;0;0;0];
Beta1 = unwrap(angle(phih))*(360/(2*pi)); % get optmal phase angle

[t,x]=ode45(@(t,x)cdx(t,x,Beta1,O),simtime,xinit);

figure('Position', [300, 300, 1200, 600]);
subplot(211)
plot(t,x(:,1)); hold on;
plot(t,x(:,3)); plot(t,x(:,5)); plot(t,x(:,7)); plot(t,x(:,9));
xlabel('Time(s)','FontSize', 15); 
ylabel('{\zeta}_{i1} {(t)}', 'FontSize', 15);
title('Coupled AHO','FontSize', 15);
axis([0 2 -1.2 1.2]); grid on;
subplot(212)
plot(t,x(:,9)); hold on;
plot(t,x(:,10));
xlabel('Time(s)','FontSize', 15); 
ylabel('{\zeta}_{51} {(t)} and {\zeta}_{52} {(t)}', 'FontSize', 15);
axis([0 2 -1.2 1.2]); grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta1 = Beta1(1); beta2 = Beta1(2); beta3 = Beta1(3); beta4 = Beta1(4); beta5 = Beta1(5);
X1 = X(1:2,:); X2 = X(3:4,:); X3 = X(5:6,:); X4 = X(7:8,:); X5 = X(9:10,:); 
U1 = U(1,:); U2 = U(2,:); U3 = U(3,:); U4 = U(4,:); U5 = U(5,:); 
h1 = [sin(beta1) cos(beta1);sin(beta1+90) cos(beta1+90)];
h2 = [sin(beta2) cos(beta2);sin(beta2+90) cos(beta2+90)];
h3 = [sin(beta3) cos(beta3);sin(beta3+90) cos(beta3+90)];
h4 = [sin(beta4) cos(beta4);sin(beta4+90) cos(beta4+90)];
h5 = [sin(beta5) cos(beta5);sin(beta5+90) cos(beta5+90)];
H = [h1;h2;h3;h4;h5];
XX = blkdiag(X1*inv(h1),X2*inv(h2),X3*inv(h3),X4*inv(h4),X5*inv(h5));
UU = blkdiag(U1*inv(h1),U2*inv(h2),U3*inv(h3),U4*inv(h4),U5*inv(h5));
YY = Cp*XX;
options=odeset('RelTol',1e-4,'Refine',5);
simtime=[0 2];
xinit=[sin(beta1);cos(beta1);sin(beta2);cos(beta2);sin(beta3);
    cos(beta3);sin(beta4);cos(beta4);sin(beta5);cos(beta5);
    0;0;0;0;0;0;0;0;0;0;0];
Beta = unwrap(angle(phih))*(360/(2*pi));
K = blkdiag(0.00009,0.00009,0.00009,0.00009,0.00009);

[tt,xx]=ode45(@(tt,xx) ddx(tt,xx,K,P,Bp,Cp,XX,UU,Beta,O,phih),simtime,xinit);

u = [];
for i = 1:size(tt)
    rho = 0.05*sin(300*tt(i));
    y = (1+rho)*Cp*xx(i,11:20)';
    u(i,:) = UU*xx(i,1:10)'- K*(y-Cp*XX*xx(i,1:10)');
end
figure('Position', [350, 300, 1200, 900]);
subplot(311)
plot(tt,(360/(2*pi))*xx(:,11)); hold on;
plot(tt,(360/(2*pi))*xx(:,13)); plot(tt,(360/(2*pi))*xx(:,15)); plot(tt,(360/(2*pi))*xx(:,17));plot(tt,(360/(2*pi))*xx(:,19));
title('Output Regulator','FontSize', 15); 
axis([0 2 -80 80]);
xlabel('Time(s)','FontSize', 15); 
ylabel('{\phi}_{i} {(t) [deg]}', 'FontSize', 15); grid on;
subplot(312)
plot(tt,u); hold on;
axis([0 2 -2*10^(-4) 2*10^(-4)]);
xlabel('Time(s)','FontSize', 15); 
ylabel('u_{i} (t) [Nm]', 'FontSize', 15); grid on;
subplot(313)
plot(tt,xx(:,21)); hold on;
xlabel('Time(s)','FontSize', 15); 
ylabel('v (t) [m/s]', 'FontSize', 15); grid on;
axis([0 2 -0.18 0.03]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (e)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xinit=[sin(beta1);cos(beta1);sin(beta2);cos(beta2);sin(beta3);
    cos(beta3);sin(beta4);cos(beta4);sin(beta5);cos(beta5);
    0;0;0;0;0;0;0;0;0;0;0];
Beta = unwrap(angle(phih))*(360/(2*pi));
K = blkdiag(0.0002,0.0002,0.0002,0.0002,0.0002);
LL = blkdiag([100;100],[100;100],[100;100],[100;100],[100;100]);
eps = 0.1;
simtime2 = [0 4];

% With disturbance without Q2
[ttt,xxx]=ode45(@(ttt,xxx) edx1(ttt,xxx,K,P,Bp,Cp,XX,UU,Beta,O,phih,LL,eps),simtime2,xinit);

% With disturbance with Q2
[tttt,xxxx]=ode45(@(tttt,xxxx) edx2(tttt,xxxx,K,P,Bp,Cp,XX,UU,Beta1,O,phih,LL,eps),simtime2,xinit);

uu = [];
d = 20;
for i = 1:size(ttt)
    rho = 0;
    y = (1+rho)*Cp*xxx(i,11:20)';
    if 2 <= ttt
        if ttt <= 2.1
            uu(i,:) = d*(UU*xxx(i,1:10)'-K*(y-Cp*XX*xxx(i,1:10)'));
        else
            uu(i,:) = UU*xxx(i,1:10)'-K*(y-Cp*XX*xxx(i,1:10)');
        end
    else
        uu(i,:) = UU*xxx(i,1:10)'-K*(y-Cp*XX*xxx(i,1:10)');    
    end
end
uuu = [];
d = 20;
for i = 1:size(tttt)
    rho = 0;
    y = (1+rho)*Cp*xxxx(i,11:20)';
    if 2 <= tttt
        if tttt <= 2.1
            uuu(i,:) = d*(UU*xxxx(i,1:10)'-K*(y-Cp*XX*xxxx(i,1:10)'));
        else
            uuu(i,:) = UU*xxxx(i,1:10)'-K*(y-Cp*XX*xxxx(i,1:10)');
        end
    else
        uuu(i,:) = UU*xxxx(i,1:10)'-K*(y-Cp*XX*xxxx(i,1:10)');    
    end
end

figure('Position', [400, 300, 1200, 900]);
subplot(311)
plot(ttt,(360/(2*pi))*xxx(:,11)); hold on;
plot(ttt,(360/(2*pi))*xxx(:,13)); plot(ttt,(360/(2*pi))*xxx(:,15)); plot(ttt,(360/(2*pi))*xxx(:,17));plot(ttt,(360/(2*pi))*xxx(:,19));
title('Output Regulator under disturbance','FontSize', 15); 
axis([0 4 -80 80]);
xlabel('Time(s)','FontSize', 15); 
ylabel('{\phi}_{i} {(t) [deg]}', 'FontSize', 15); grid on;
subplot(312)
plot(ttt,uu); hold on;
axis([0 4 -2*10^(-4) 2*10^(-4)]);
xlabel('Time(s)','FontSize', 15); 
ylabel('u_{i} (t) [Nm]', 'FontSize', 15); grid on;
subplot(313)
plot(ttt,xxx(:,21)); hold on;
xlabel('Time(s)','FontSize', 15); 
ylabel('v (t) [m/s]', 'FontSize', 15); grid on;
axis([0 4 -0.4 0.9]);

figure('Position', [450, 300, 1200, 900]);
subplot(311)
plot(tttt,(360/(2*pi))*xxxx(:,11)); hold on;
plot(tttt,(360/(2*pi))*xxxx(:,13)); plot(tttt,(360/(2*pi))*xxxx(:,15)); plot(tttt,(360/(2*pi))*xxxx(:,17));plot(tttt,(360/(2*pi))*xxxx(:,19));
title('Output Regulator with Feedback','FontSize', 15); 
axis([0 4 -80 80]);
xlabel('Time(s)','FontSize', 15); 
ylabel('{\phi}_{i} {(t) [deg]}', 'FontSize', 15); grid on;
subplot(312)
plot(tttt,uuu); hold on;
axis([0 4 -2*10^(-4) 2*10^(-4)]);
xlabel('Time(s)','FontSize', 15); 
ylabel('u_{i} (t) [Nm]', 'FontSize', 15); grid on;
subplot(313)
plot(tttt,xxxx(:,21)); hold on;
xlabel('Time(s)','FontSize', 15); 
ylabel('v (t) [m/s]', 'FontSize', 15); grid on;
axis([0 4 -0.4 0.9]);





