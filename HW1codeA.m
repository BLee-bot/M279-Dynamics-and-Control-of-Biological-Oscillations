% HW1 (a) BANSEOK LEE
clear all
options=odeset('RelTol',1e-4,'Refine',5);
simtime=[0 4000];
xinit=[0;0;0;0;0;0;0;0;0];
[t,x]=ode45(@ahf,simtime,xinit);
[t2,x2]=ode45(@ahf2,simtime,xinit);
[t3,x3]=ode45(@ahf3,simtime,xinit);
[t4,x4]=ode45(@ahf4,simtime,xinit);
[t5,x5]=ode45(@ahf5,simtime,xinit);
[t6,x6]=ode45(@ahf6,simtime,xinit);

figure(1)
subplot(411)
plot(t/1000,x(:,1)); 
hold on;
plot(t/1000,x(:,4)); 
plot(t/1000,x(:,7)); 
title('sigma = 1.50'); xlabel('Time(s)'); ylabel('v_i');
subplot(412)
plot(t2/1000,x2(:,1)); 
hold on;
plot(t2/1000,x2(:,4)); 
plot(t2/1000,x2(:,7)); 
title('sigma = 1.89'); xlabel('Time(s)'); ylabel('v_i');
subplot(413)
plot(t3/1000,x3(:,1)); 
hold on;
plot(t3/1000,x3(:,4)); 
plot(t3/1000,x3(:,7)); 
title('sigma = 1.90'); xlabel('Time(s)'); ylabel('v_i');
subplot(414)
plot(t4/1000,x4(:,1)); 
hold on;
plot(t4/1000,x4(:,4)); 
plot(t4/1000,x4(:,7)); 
title('sigma = 2.00'); xlabel('Time(s)'); ylabel('v_i');



figure(2)
subplot(211)
plot(t5/1000,x5(:,1)); 
hold on;
plot(t5/1000,x5(:,4)); 
plot(t5/1000,x5(:,7)); 
title('sigma = 2.00 / +3 Impulse applied u1 at 3.3sec')
xlabel('Time(s)'); ylabel('v_i');

subplot(212)
plot(t6/1000,x6(:,1)); 
hold on;
plot(t6/1000,x6(:,4)); 
plot(t6/1000,x6(:,7)); 
title('sigma = 2.00 / +0.3 Impulse applied u1 at 3.3sec')
xlabel('Time(s)'); ylabel('v_i');








