clear; clc; close all;
%% Q1.1 Kinematics
syms r_m theta_l theta_r dot_theta_l dot_theta_r ddot_theta_l ddot_theta_r

x = r_m/2*(theta_l+theta_r);
y = r_m/2*(theta_l-theta_r);
theta_matrix = [theta_l; theta_r];
theta_dot_matrix = [dot_theta_l; dot_theta_r];
theta_ddot_matrix = [ddot_theta_l; ddot_theta_r];
xy_matrix = [x;y]
xy_dot_matrix = subs(xy_matrix, [theta_l, theta_r], [dot_theta_l, dot_theta_r])
xy_ddot_matrix = subs(xy_matrix, [theta_l, theta_r], [ddot_theta_l, ddot_theta_r])

H = r_m/2 * [1 1; 1 -1];
%% Q1.2 Dynamics
syms J_m m mu_l_x mu_l_y d_roll i_l i_r K_m u_l u_r dot_i_l dot_i_r R_m L_m

J_eq = J_m + 1/2*m*r_m^2;
d_plus = r_m^2/4*(mu_l_x+mu_l_y) + d_roll;
d_minus = r_m^2/4*(mu_l_x-mu_l_y);

i_matrix = [i_l; i_r];
i_dot_matrix = [dot_i_l; dot_i_r];
u_matrix = [u_l; u_r];

J_eq*theta_ddot_matrix + [d_plus, d_minus; d_minus, d_plus]*theta_dot_matrix == K_m*i_matrix

u_matrix == R_m*i_matrix + L_m*i_dot_matrix + K_m*theta_dot_matrix

A = sym(zeros(6,6));
A(1,3) = 1;
A(2,4) = 1;
A(3,:) = [0, 0, -d_plus/J_eq, -d_minus/J_eq, K_m/J_eq, 0];
A(4,:) = [0, 0, -d_minus/J_eq, -d_plus/J_eq, 0, K_m/J_eq];
A(5,:) = [0, 0, -K_m/L_m, 0, -R_m/L_m, 0];
A(6,:) = [0, 0, 0, -K_m/L_m, 0, -R_m/L_m];
B = sym(zeros(6,2));
B(5,1) = 1/L_m;
B(6,2) = 1/L_m;

X = [theta_l; theta_r; dot_theta_l; dot_theta_r; i_l; i_r]
X_dot = [dot_theta_l; dot_theta_r; ddot_theta_l; ddot_theta_r; dot_i_l; dot_i_r]
U = u_matrix

X_dot == A*X + B*U
%% Q1.3 Decoupled dynamics
i_x = i_l + i_r;
i_y = i_l - i_r;
u_x = u_l + u_r;
u_y = u_l - u_r;
Z = [xy_matrix(1); xy_dot_matrix(1); i_x; xy_matrix(2); xy_dot_matrix(2); i_y]
V = [u_x; u_y]
T = [r_m/2 r_m/2 0 0 0 0;
    0 0 r_m/2 r_m/2 0 0;
    0 0 0 0 1 1;
    r_m/2 -r_m/2 0 0 0 0;
    0 0 r_m/2 -r_m/2 0 0;
    0 0 0 0 1 -1];
simplify(Z - T*X)
S = [1 1; 1 -1];

A_decoupled = T*A/T
B_decoupled = T*B/S

A_x = A_decoupled(1:3,1:3)
A_y = A_decoupled(4:6,4:6)
B_x = B_decoupled(1:3,1)
B_y = B_decoupled(4:6,2)
%% Q1.4 Drawing Ellipse
syms t a b omega x_c y_c
x_d = x_c + a*cos(omega*t);
y_d = y_c + b*sin(omega*t);
xy_d_matrix = [x_d; y_d]
xy_d_dot_matrix = [jacobian(x_d,t); jacobian(y_d,t)]

alpha = r_m*K_m/(2*J_eq);
beta_x= (d_plus+d_minus)/J_eq;
beta_y= (d_plus-d_minus)/J_eq;

i_x_d = (jacobian(xy_d_dot_matrix(1),t) + beta_x*xy_d_dot_matrix(1))/alpha;
i_y_d = (jacobian(xy_d_dot_matrix(2),t) + beta_y*xy_d_dot_matrix(2))/alpha;
i_d_matrix = [i_x_d; i_y_d]

u_xy_d_matrix = L_m*jacobian(i_d_matrix, t) + R_m*i_d_matrix + 2*K_m/r_m*xy_d_dot_matrix
u_lr_d_matrix = 1/2*S*u_xy_d_matrix
% These should be our inputs; now to test the outcome from the initial conditions
x_d_0 = subs(x_d, t, 0);
x_d_dot = jacobian(x_d, t);
x_d_dot_0 = subs(x_d_dot, t, 0);
i_x_d_0 = subs(i_x_d, t, 0);
y_d_0 = subs(y_d, t, 0);
y_d_dot = jacobian(y_d, t);
y_d_dot_0 = subs(y_d_dot, t, 0);
i_y_d_0 = subs(i_y_d, t, 0);
Z_0 = [x_d_0; x_d_dot_0; i_x_d_0; y_d_0; y_d_dot_0; i_y_d_0]

X_0 = inv(T)*Z_0

parameters = [r_m; J_m; m; mu_l_x; mu_l_y; d_roll; K_m; R_m; L_m; a; b; omega; x_c; y_c];
parameters_values = [0.010; 1.0e-4; 0.50; 0.50; 0.70; 1.0e-5; 0.050; 2.0; 0.010; 0.050; 0.035; 1.0; 0.15; 0.15];
A_sim = double(subs(A, parameters, parameters_values));
B_sim = double(subs(B, parameters, parameters_values));
X_0_sim = double(subs(X_0, parameters, parameters_values));
U_d_sim = subs(u_lr_d_matrix, parameters, parameters_values);

n_steps = 5001;
t_sim = linspace(0, 2*pi, n_steps);
X_sim = zeros(size(X,1), n_steps);
X_sim(:,1) = X_0_sim;
for i = 1:n_steps - 1
    U_sim = double(subs(U_d_sim, t, t_sim(i)));
    X_dot_sim = A_sim*X_sim(:,i) + B_sim*U_sim;
    X_sim(:, i+1) = X_sim(:, i) + X_dot_sim*(t_sim(i+1)-t_sim(i));
end

x_sim = double(subs(r_m/2*(X_sim(1,:)+X_sim(2,:)), r_m, 0.010));
y_sim = double(subs(r_m/2*(X_sim(1,:)-X_sim(2,:)), r_m, 0.010));

x_d_sim = subs(x_d, [x_c, a, omega], [0.15, 0.050, 1.0]);
x_d_sim = subs(x_d_sim, t, linspace(0, 2*pi, 101));
y_d_sim = subs(y_d, [y_c, b, omega], [0.15, 0.035, 1.0]);
y_d_sim = subs(y_d_sim, t, linspace(0, 2*pi, 101));

figure;
plot(x_sim, y_sim, 'g-', 'LineWidth', 1.5);
hold on
plot(x_d_sim, y_d_sim, 'k--', 'LineWidth', 1.5);
hold off
grid on;
xlim([0.075 0.225])
ylim([0.075 0.225])
xlabel('x [m]');
ylabel('y [m]');
title('Q1.4: Desired and simulated ellipse trajectory');
legend('Simulated trajectory', 'Desired trajectory', 'Location', 'best');

for i = 1:n_steps
    U_sim(:,i) = double(subs(U_d_sim, t, t_sim(i)));
end

figure;
plot(t_sim, U_sim(1,:), 'LineWidth', 1.5);
hold on;
plot(t_sim, U_sim(2,:), 'LineWidth', 1.5);
grid on;
xlim([0 2*pi])
xlabel('Time [s]');
ylabel('Motor voltage [V]');
title('Q1.4: Required motor voltages');
legend('Left motor', 'Right motor', 'Location', 'best');
%% Q1.4 Drawing infinity sign
syms t a b omega x_c y_c
x_d = x_c + a*cos(omega*t);
y_d = y_c + b*sin(omega*t)*cos(omega*t);
xy_d_matrix = [x_d; y_d]
xy_d_dot_matrix = [jacobian(x_d,t); jacobian(y_d,t)]

alpha = r_m*K_m/(2*J_eq);
beta_x= (d_plus+d_minus)/J_eq;
beta_y= (d_plus-d_minus)/J_eq;

i_x_d = (jacobian(xy_d_dot_matrix(1),t) + beta_x*xy_d_dot_matrix(1))/alpha;
i_y_d = (jacobian(xy_d_dot_matrix(2),t) + beta_y*xy_d_dot_matrix(2))/alpha;
i_d_matrix = [i_x_d; i_y_d]

u_xy_d_matrix = L_m*jacobian(i_d_matrix, t) + R_m*i_d_matrix + 2*K_m/r_m*xy_d_dot_matrix
u_lr_d_matrix = 1/2*S*u_xy_d_matrix
% These should be our inputs; now to test the outcome from the initial conditions
x_d_0 = subs(x_d, t, 0);
x_d_dot = jacobian(x_d, t);
x_d_dot_0 = subs(x_d_dot, t, 0);
i_x_d_0 = subs(i_x_d, t, 0);
y_d_0 = subs(y_d, t, 0);
y_d_dot = jacobian(y_d, t);
y_d_dot_0 = subs(y_d_dot, t, 0);
i_y_d_0 = subs(i_y_d, t, 0);
Z_0 = [x_d_0; x_d_dot_0; i_x_d_0; y_d_0; y_d_dot_0; i_y_d_0]

X_0 = inv(T)*Z_0

parameters = [r_m; J_m; m; mu_l_x; mu_l_y; d_roll; K_m; R_m; L_m; a; b; omega; x_c; y_c];
parameters_values = [0.010; 1.0e-4; 0.50; 0.50; 0.70; 1.0e-5; 0.050; 2.0; 0.010; 0.050; 0.035; 1.0; 0.15; 0.15];
A_sim = double(subs(A, parameters, parameters_values));
B_sim = double(subs(B, parameters, parameters_values));
X_0_sim = double(subs(X_0, parameters, parameters_values));
U_d_sim = subs(u_lr_d_matrix, parameters, parameters_values);

n_steps = 5001;
t_sim = linspace(0, 2*pi, n_steps);
X_sim = zeros(size(X,1), n_steps);
X_sim(:,1) = X_0_sim;
for i = 1:n_steps - 1
    U_sim = double(subs(U_d_sim, t, t_sim(i)));
    X_dot_sim = A_sim*X_sim(:,i) + B_sim*U_sim;
    X_sim(:, i+1) = X_sim(:, i) + X_dot_sim*(t_sim(i+1)-t_sim(i));
end

x_sim = double(subs(r_m/2*(X_sim(1,:)+X_sim(2,:)), r_m, 0.010));
y_sim = double(subs(r_m/2*(X_sim(1,:)-X_sim(2,:)), r_m, 0.010));

x_d_sim = subs(x_d, [x_c, a, omega], [0.15, 0.050, 1.0]);
x_d_sim = subs(x_d_sim, t, linspace(0, 2*pi, 101));
y_d_sim = subs(y_d, [y_c, b, omega], [0.15, 0.035, 1.0]);
y_d_sim = subs(y_d_sim, t, linspace(0, 2*pi, 101));

figure;
plot(x_sim, y_sim, 'g-', 'LineWidth', 1.5);
hold on
plot(x_d_sim, y_d_sim, 'k--', 'LineWidth', 1.5);
hold off
grid on;
xlim([0.075 0.225])
ylim([0.075 0.225])
xlabel('x [m]');
ylabel('y [m]');
title('Q1.4: Desired and simulated infinity sign trajectory');
legend('Simulated trajectory', 'Desired trajectory', 'Location', 'best');

for i = 1:n_steps
    U_sim(:,i) = double(subs(U_d_sim, t, t_sim(i)));
end

figure;
plot(t_sim, U_sim(1,:), 'LineWidth', 1.5);
hold on;
plot(t_sim, U_sim(2,:), 'LineWidth', 1.5);
grid on;
xlim([0 2*pi])
xlabel('Time [s]');
ylabel('Motor voltage [V]');
title('Q1.4: Required motor voltages');
legend('Left motor', 'Right motor', 'Location', 'best');