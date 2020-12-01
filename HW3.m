% HW3 BANSEOK LEE 605351891
clear all
close all

% Define Constants
n = 5; m = 0.001; l = 0.1; v0 = -0.15;
cni = 0.009; cti = 0.0006; ct0 = 0.0006;
li = l/(2*n+2); mi = m/(n+1); l0 = l/(2*n+2); m0 = m/(n+1);

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
S = (Lamda - Lamda')/2;

% History array
lamda_hist = [];
w_hist = [];
u_hat_hist = [];

% Line search for w
k=1;
for w = 1:0.1:100
    P_w = inv(-J*w^2+D*(1i*w)+v0*Lamda)*B;
    X_w = (1/4)*(1i*w)*(B'*P_w-P_w'*B);
    Y_w = (1/(2*(ct0+e'*Ct*e)*v0))*P_w'*(1i*w*S - v0*C0)*P_w;
    [V,E] = eig(X_w, Y_w); % Generalized Eigenvalues, Eigenvectors
    R_idx = abs(imag(diag(E)))<0.0000000000001; % Threshold for real number
    R = E*R_idx;
    idx = [];
    for l = 1:n
        if R_idx(l) == 0
            idx(l) = 0;
        else
            if eig(X_w-R(l)*Y_w)>=0
                idx(l) = 1;
            else
                idx(l) = 0;
            end
        end
    end
    [Max, argmax] = max(E*idx');
    v = V(:,argmax); 
    lamda_hist(k) = (real(Max)>0)*real(Max);
    w_hist(k) = w;
    u_hat_hist(:,k) = (1/sqrt(v'*Y_w*v))*v;
    k=k+1;
end

% Delete all noise
w_hist(lamda_hist<0.0001) = [];
u_hat_hist(:,lamda_hist<0.0001) = [];
lamda_hist(lamda_hist<0.0001) = [];

% Find optimal input and frequency
[~, argmin] = min(lamda_hist);
opt_w = w_hist(argmin);
u_hat = u_hat_hist(:,argmin);

% Normalize u hat, theta hat
u_hat = u_hat*exp(1j*0.3499*2*pi);
Pwo = inv(-J*opt_w^2+D*(1i*opt_w)+v0*Lamda)*B;
theta_hat = Pwo*u_hat;

% Show the normalized input signal, theta,and gait through time 
t = 0:0.01:1;
u_t = real(u_hat*exp(1j*(opt_w*t)));
theta_t = real(theta_hat*exp(1j*(opt_w*t)));
gait = B'*theta_t;

% Figure 1 Power optimal gait
figure('Position', [300, 300, 1200, 800]);
subplot(2,2,1)
plot(t,1000*1000*u_t);
axis([0 1 -150 150]);
grid on;
ylabel('u(t) [mN.mm]', 'FontSize', 20);
xlabel('Time [s]', 'FontSize', 20);
subplot(2,2,2)
plot(t,(360/(2*pi))*gait);
axis([0 1 -50 50]);
grid on;
ylabel('{\phi}(t) [deg]', 'FontSize', 20);
xlabel('Time [s]', 'FontSize', 20);
gaithat = B'*theta_hat;
idx = 1:n;
subplot(2,2,3)
plot(idx, (360/(2*pi))*angle(u_hat)); hold on;
plot(idx, (360/(2*pi))*unwrap(angle(gaithat)));
axis([0 6 -360 100]);
grid on;
ylabel('{\angle} {u}_{i} (b), {\angle} {\phi}_{i} (r) [deg]', 'FontSize', 20);
xlabel('Link Index', 'FontSize', 20);
subplot(2,2,4)
plot(idx, 1000*500*abs(u_hat)); hold on;
plot(idx, (360/(2*pi))*abs(gaithat));
axis([0 6 0 70]);
grid on;
ylabel('|ui|/2 (b) [mN.mm], |fi| (r) [deg]', 'FontSize', 20);
xlabel('Link Index', 'FontSize', 20);

% Figure 2 Optimal frequency f_opt = 2.51 Hz, optimal power 0.253 mW.
figure('Position', [300, 300, 600, 600]);
x = w_hist./(2*pi); y = lamda_hist*1000;
semilogx(x,y); hold on;
[~,imn] = min(y);
semilogx(x(imn),y(imn),'pr')
ylabel('Input Power [mW]', 'FontSize', 20);
xlabel('{\omega} [Hz]', 'FontSize', 20);
grid on;
axis([0 15 0 1]);
x(imn)
y(imn)
% Figure 3 Snapshots of the body under the optimal gait
figure('Position', [300, 300, 1200, 200]);
hold on;
for T = (0.1:0.1:1)*2*pi/opt_w
    G = B'*imag(theta_hat*exp(1j*(opt_w*T)));
    Locomotion = [];
    Locomotion(1,:) = [-0.15*T-2*l 0];
    Locomotion(2,:) = [-0.15*T 0];
    Locomotion(3,:) = Locomotion(2,:) + 2*li*[cos(G(1)) sin(G(1))];
    Locomotion(4,:) = Locomotion(3,:) + 2*li*[cos(G(2)) sin(G(2))];
    Locomotion(5,:) = Locomotion(4,:) + 2*li*[cos(G(3)) sin(G(3))];
    Locomotion(6,:) = Locomotion(5,:) + 2*li*[cos(G(4)) sin(G(4))];
    Locomotion(7,:) = Locomotion(6,:) + 2*li*[cos(G(5)) sin(G(5))];
    plot(Locomotion(:,1),Locomotion(:,2),'b')
end
axis([-0.07 0.07 -0.015 0.015]);
(360/(2*pi))*unwrap(angle(gaithat))
