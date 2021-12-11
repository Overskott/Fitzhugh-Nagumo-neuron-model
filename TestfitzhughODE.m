% Modified from TestConvergenceODE.m from the Oslomet ACIT4310 course
% Compares the fitznag.m versus Euler, Heuns and RK4

% Values of h: 1, 1/2, 1/2^2, ..., 1/2^10.
h=2.^[0:-1:-6];
Nh=length(h);
tmax=10;
EE=zeros(Nh,1);
EH=zeros(Nh,1);
ERK=zeros(Nh,1);

IC=[0 0];

for k=1:Nh
    t=[0:h(k):tmax];
    % Exact solution:
    [time, YExact]=ode45(@fitznag,t,IC);

    Y = (Euler(t, IC, @fitznag))';
    EE(k)=mean(abs(Y(:,1)- YExact(:,1)));
    
    Y = (Heun(t, IC, @fitznag))';
    EH(k)=mean(abs(Y(:,1)- YExact(:,1)));
    
    Y = (RungeKutta4(t, IC, @fitznag))';
    ERK(k)=mean(abs(Y(:,1)- YExact(:,1)));
end

hold on
plot(log(h),log(EE),'-o','linewidth',1)
te1=['Euler'];

plot(log(h),log(EH),'-o','linewidth',1)
te3=['Heun'];

plot(log(h),log(ERK),'-o','linewidth',1)
te5=['Runge-Kutta4'];


xlabel('ln(h)')
ylabel('ln(|error|)')
set(gca,'fontsize',12)
set(legend,'location', 'best')

legend(te1,te3,te5)
hold off




