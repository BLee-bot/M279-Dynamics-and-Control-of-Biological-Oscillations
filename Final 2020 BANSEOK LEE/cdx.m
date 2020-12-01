
function cdx = cdx(t,x,Beta, O)
    
    x1 = x(1:2,:); x2 = x(3:4,:); x3 = x(5:6,:); x4 = x(7:8,:); x5 = x(9:10,:);
    beta1 = Beta(1); beta2 = Beta(2); beta3 = Beta(3); beta4 = Beta(4); beta5 = Beta(5);
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
    cdx = (blkdiag(A1,A2,A3,A4,A5) + R)*x;
end
