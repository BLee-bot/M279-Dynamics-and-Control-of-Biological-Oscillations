% HW2 (a) BANSEOK LEE 
clear all
options=odeset('RelTol',1e-4,'Refine',5);
simtime=[0 10];
xinit=[3;0;1;0];
[t,x]=ode45(@dAHOa,simtime,xinit);

figure(1)
plot(t,x(:,1)); hold on;
plot(t,x(:,2));
plot(t,x(:,3)); 
plot(t,x(:,4));
title('x1, x2 over t'); xlabel('Time(s)'); ylabel('x1,x2');

figure(2)
subplot(121)
plot(x(:,1),x(:,2)); axis([-3 3  -3 3]); axis square;
title('Plot x1 = (y1, z1)'); xlabel('y1'); ylabel('z1');
subplot(122)
plot(x(:,3),x(:,4)); axis([-3 3  -3 3]); axis square;
title('Plot x2 = (y2, z2)'); xlabel('y2'); ylabel('z2');






