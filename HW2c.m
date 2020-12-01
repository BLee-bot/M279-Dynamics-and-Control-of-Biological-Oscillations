% HW2 (c) BANSEOK LEE 
clear all
options=odeset('RelTol',1e-4,'Refine',5);
simtime=[0 30];
simtime2 = 0:0.03:30;
xinit=[3;0;1;0];
XX = []; EE = [];
i = 0; k=0;

% Calculate ode for various (delta1, delta2) pairs
for delta1 = flip(-3:0.5:3)
    j = 0;
    i = i+1;
    for delta2 = -3:0.5:3
        j = j+1;
        [t,x]=ode45(@(t,x)dAHOc(t,x,delta1,delta2),simtime2,xinit);
        XX(:,:,i,j)=x;
    end
end

% Plot
figure(1)
title('Synchronized Error graph for (delta1, delta2)');
for n = 1:13
    for m = 1:13
        k=k+1;
        subplot(13,13,k)        
        plot(t,XX(:,1,n,m)-XX(:,3,n,m)); hold on;
        plot(t,XX(:,2,n,m)-XX(:,4,n,m)); hold off;
        set(gca,'xtick',[],'ytick',[])
        EE(1,n,m) = 0; EE(2,n,m) = 0;
        for l = 800:900
            EE(1,n,m) = EE(1,n,m)+((XX(l,1,n,m)-XX(l,3,n,m))^2)*(t(l+1)-t(l));
            EE(2,n,m) = EE(2,n,m)+((XX(l,2,n,m)-XX(l,4,n,m))^2)*(t(l+1)-t(l));
        end
        title(mean([EE(1,n,m),EE(2,n,m)]));
    end
end

