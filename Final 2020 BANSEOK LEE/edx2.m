
function edx2 = edx2(tttt,xxxx,K,P,Bp,Cp,XX,UU,Beta1,O,phih,LL,eps)
    
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
    
    edx2 = [];
    x1 = xxxx(1:2,:);  % zeta1
    x2 = xxxx(3:4,:);  % zeta2
    x3 = xxxx(5:6,:);  % zeta3
    x4 = xxxx(7:8,:);  % zeta4
    x5 = xxxx(9:10,:); % zeta5
    x0 = xxxx(1:10,:); % Zeta
    x6 = xxxx(11:20,:);% state x
    v = xxxx(21,:);    % v 
    
    % AHO
    beta1 = Beta1(1); beta2 = Beta1(2); beta3 = Beta1(3); beta4 = Beta1(4); beta5 = Beta1(5);
    mu = 200;
    A1 = -O + eye(2)*(mu)*(1-norm(x1)^2);
    A2 = -O + eye(2)*(mu)*(1-norm(x2)^2);
    A3 = -O + eye(2)*(mu)*(1-norm(x3)^2);
    A4 = -O + eye(2)*(mu)*(1-norm(x4)^2);
    A5 = -O + eye(2)*(mu)*(1-norm(x5)^2);
    R1 = [cos(beta1) -sin(beta1);sin(beta1) cos(beta1)];
    R2 = [cos(beta2) -sin(beta2);sin(beta2) cos(beta2)];
    R3 = [cos(beta3) -sin(beta3);sin(beta3) cos(beta3)];
    R4 = [cos(beta4) -sin(beta4);sin(beta4) cos(beta4)];
    R5 = [cos(beta5) -sin(beta5);sin(beta5) cos(beta5)];
    R = [-eye(2) R1'*R2 zeros(2) zeros(2) zeros(2);
        R2'*R1 -2*eye(2) R2'*R3 zeros(2) zeros(2);
        zeros(2) R3'*R2 -2*eye(2) R3'*R4 zeros(2);
        zeros(2) zeros(2) R4'*R3 -2*eye(2) R4'*R5;
        zeros(2) zeros(2) zeros(2) R5'*R4 -eye(2)];
    
    % get feedback u
    w = O(2,1);
    h1 = [sin(beta1) cos(beta1);sin(beta1+90) cos(beta1+90)];
    h2 = [sin(beta2) cos(beta2);sin(beta2+90) cos(beta2+90)];
    h3 = [sin(beta3) cos(beta3);sin(beta3+90) cos(beta3+90)];
    h4 = [sin(beta4) cos(beta4);sin(beta4+90) cos(beta4+90)];
    h5 = [sin(beta5) cos(beta5);sin(beta5+90) cos(beta5+90)];
    uhat = inv(B)*(-J*w^2+1i*w*D+v*Lamda)*inv(B')*phih;
    Gamma = 0.5*[-1i 1; 1i 1];
    U = [uhat conj(uhat)]*Gamma;
    U1 = U(1,:); U2 = U(2,:); U3 = U(3,:); U4 = U(4,:); U5 = U(5,:); 
    UU = blkdiag(U1*inv(h1),U2*inv(h2),U3*inv(h3),U4*inv(h4),U5*inv(h5));
    rho = 0; % zero noise
    y = (1+rho)*Cp*x6;
    
    % add disturbance
    d = 20;
    if 2 <= tttt 
        if tttt <= 2.1
            u = d*(UU*x0-K*(y-Cp*XX*x0));
        else
            u = UU*x0-K*(y-Cp*XX*x0);
        end
    else
        u = UU*x0-K*(y-Cp*XX*x0);
    end
    
    % get Zeta dot with feedback L
    edx2(1:10,:) = (blkdiag(A1,A2,A3,A4,A5) + R)*x0 + eps*LL*(y-Cp*XX*x0); 
    
    % get state dot from nonlinear plant dynamics
    AA = [zeros(n) eye(n);-B'*inv(J)*v*Lamda*inv(B') -B'*inv(J)*D*inv(B')];
    AAA = P*AA*inv(P);
    edx2(11:20,:) = AAA*x6 +Bp*u;
    
    % get v dot 
    thetadot = inv(B')*[xxxx(12,:);xxxx(14,:);xxxx(16,:);xxxx(18,:);xxxx(20,:)];
    theta = inv(B')*[xxxx(11,:);xxxx(13,:);xxxx(15,:);xxxx(17,:);xxxx(19,:)];
    dtheta = ct0 +e'*Ct*e+theta'*C0*theta;
    edx2(21,:) = (-dtheta*v-thetadot'*Lamda*theta)/m;
   
end
