% HW1 (c) BANSEOK LEE
clear all
options=odeset('RelTol',1e-4,'Refine',5);
simtime=[0 40000];
xinit=[0;0;0];
[t,x]=ode45(@chf2,simtime,xinit);
[t2,x2]=ode45(@chf3,simtime,xinit);
[t3,x3]=ode45(@chf4,simtime,xinit);
[t4,x4]=ode45(@chf5,simtime,xinit);
[t5,x5]=ode45(@chf6,simtime,xinit);
[t6,x6]=ode45(@chf7,simtime,xinit);

figure(1)
subplot(411)
plot(t/1000,x(:,1)); 
hold on;
plot(t/1000,x(:,2)); 
plot(t/1000,x(:,3)); 
xlabel('Time(s)'); ylabel('q_i'); title('sigma = 1.5');
axis([0 20 -1 1])
subplot(412)
plot(t2/1000,x2(:,1)); 
hold on;
plot(t2/1000,x2(:,2)); 
plot(t2/1000,x2(:,3)); 
xlabel('Time(s)'); ylabel('q_i'); title('sigma = 1.90');
axis([0 20 -1 1])
subplot(413)
plot(t3/1000,x3(:,1)); 
hold on;
plot(t3/1000,x3(:,2)); 
plot(t3/1000,x3(:,3)); 
xlabel('Time(s)'); ylabel('q_i'); title('sigma = 2.10');
axis([0 20 -1 1])
subplot(414)
plot(t4/1000,x4(:,1)); 
hold on;
plot(t4/1000,x4(:,2)); 
plot(t4/1000,x4(:,3)); 
xlabel('Time(s)'); ylabel('q_i'); title('sigma = 2.30');
axis([0 20 -1 1])

figure(2)
subplot(211)
plot(t5/1000,x5(:,1)); 
hold on;
plot(t5/1000,x5(:,2)); 
plot(t5/1000,x5(:,3)); 
xlabel('Time(s)'); ylabel('q_i'); title('sigma = 2.5');
axis([0 4 -2 2])
subplot(212)
plot(t6/1000,x6(:,1)); 
hold on;
plot(t6/1000,x6(:,2)); 
plot(t6/1000,x6(:,3)); 
xlabel('Time(s)'); ylabel('q_i'); title('sigma = 7.0');
axis([0 4 -5 5])











