% Script to solve, plot and analyze the Fitzhugh-Nagumo model.
% Created 2021, Sebastian T. Overskott (s331402@oslomet.no)

close all % Close any open windows

% Parameters
% From https://en.wikipedia.org/w/index.php?title=FitzHugh%E2%80%93Nagumo_model&oldid=1037464531
R = 1;
I_app = 0.5;
a = 0.7;
b = 0.8;
tau = 12.5;

% Initial conditions
v_0 = 0;
w_0 = 0;

IC = [v_0, w_0];

% Numerics
t = 0:0.01:200;
n = length(t);
v = -2.5:0.01:2.5;


% ------- ODE45 Solution -------------------------------%
[time, VW] = fitznagODE(R, I_app, a, b, tau, IC, t);

v_result=VW(:,1);
w_result=VW(:,2);

% % plot
f1=figure;
figure(f1);
hold on
title("Voltage for V and W")
xlabel("Time")
ylabel("Voltage")
grid on
xline(0)
yline(0)
p_1(1)= plot(t, v_result);
p_1(2)= plot(t, w_result);
%legend(p_1,'V', 'W');
hold off

% ------- Stabel states ---------------------------%

% Coefficient vector for v'
v_coef = [-1/3 0  (b-1)/b -(a/b)+R*I_app];
v_roots = roots(v_coef);

v_real_bool = v_roots == real(v_roots); % choose only real elements
v_real = v_roots(v_real_bool); % Extract real elements

w_roots=zeros(length(v_real), 1);

% Coefficient vector for w'
for i = 1:length(v_real)
    w_coef = [-b a+v_real(i)];
    w_roots(i) = roots(w_coef);

    w_real_bool = w_roots == real(w_roots); % choose only real elements
    w_real = w_roots(w_real_bool); % Extract real elements
end

% Jacobian values

for i = 1:length(v_real)
    f_v = 1-(v_real(i)^2);
    f_w = -1;
    g_v = 1/tau;
    g_w = -b/tau;

    % Jacobian matrix
    J = [f_v f_w; g_v g_w];

    %Jacobian eigenvalues
    J_eig = eig(J);
end

% ------- Nullclines, SS and Phase line plot -----------%


nullcline_v = @(v) v-(v.^3)/3+R*I_app;
nullcline_w = @(v)(v+a)/b;

v_nc = nullcline_v(v);
w_nc = nullcline_w(v);

% Plot
f2=figure;
figure(f2);
hold on
grid on
title("Nullclines")
xlabel('V')
ylabel('W')
xline(0)
yline(0)


p_2(1)= plot(v, v_nc);
p_2(2)= plot(v, w_nc);
p_2(3)= plot(v_real, w_real, 'rx');
p_2(4)=plot(v_result,w_result);
legend(p_2,"Nullcline V'", "Nullcline W'","SS","Phase plot");
axis([-2.5 2.5 -1 3])
hold off



% ------- Bifurcation ---------------%

%I_app range
I = 0:0.001:2;

index = 1;

J_eig = length(I);
eigenvalues_1 = length(I);
eigenvalues_2 = length(I);
J_trace = length(I);
J_det = length(I);

for i = I
    
    v_coef_I = [-1/3 0 (b-1)/b -(a/b)+R*i];
    v_roots_I = roots(v_coef_I);

    v_real_bool_I = v_roots_I == real(v_roots_I); % choose only real elements
    v_real_I = v_roots_I(v_real_bool_I); % Extract real elements

    
    % Jacobian values
    f_v = 1-(v_real_I^2);
    f_w = -1;
    g_v = 1/tau;    
    g_w = -b/tau;

    % Jacobian matrix
    J = [f_v f_w; g_v g_w];

    %Jacobian eigenvalues
    J_eig = eig(J);

    % trace
    J_trace(index) = trace(J);

    % determinant
    J_det(index) = det(J);
    
    eigenvalues_1(index) = real(J_eig(1));
    eigenvalues_2(index) = real(J_eig(2));

    index = index + 1;
end

%Finding limit for I_app and stability
TraceMinimaIndex = islocalmin(abs(J_trace)); % Vector with local minima
I_value_indexed = TraceMinimaIndex.*I; % Extract corresponding I value
I_value_indexed( I_value_indexed == 0 ) = []; % Remove zero-elements
I_min = zeros(length(I_value_indexed), 1);

% plot
f3= figure;
figure(f3)
hold on
grid on
xline(0)
yline(0)
title("Bifurcation")
xlabel("I_{app}")
p_3(1)=plot(I, eigenvalues_1);
p_3(2)=plot(I, eigenvalues_2);
p_3(3)=plot(I, J_trace);
p_3(4)=plot(I, J_det);
p_3(5)=plot(I_value_indexed, I_min, 'rx');
legend(p_3, "Eig_1", "Eig_2", "Trace", "Determinant", "Change of SS")


% ------ Numerical solutions -------- %

Euler_result = Euler(t, IC, @fitznag);
Heun_result = Heun(t, IC, @fitznag);
RungeKutta4_result = RungeKutta4(t, IC, @fitznag);

% % Plot
f4 = figure;
figure(f4)
grid on
hold on
title("Comparison of different schemes with n="+num2str(n))
xlabel('time')
ylabel('V')
p_4(1)=plot(t, Euler_result(1,:), 'LineWidth', 1);
p_4(2)=plot(t, Heun_result(1,:) , 'LineWidth', 1);
p_4(3)=plot(t, RungeKutta4_result(1,:), 'LineWidth', 1);
% axis([0 10 0 inf])
legend(p_4, 'Euler', 'Heun', 'RK4');
hold off 

f5=figure;
figure(f5)
hold on
title("Comparison of schemes")
TestfitzhughODE
hold off