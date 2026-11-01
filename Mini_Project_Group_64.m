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

x_sim = double(subs(r_m/2*(X_sim(1,:) + X_sim(2,:)), r_m, 0.010));
y_sim = double(subs(r_m/2*(X_sim(1,:) - X_sim(2,:)), r_m, 0.010));

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
%% Q1.5 Model nozzle
syms T_n i_n v_n dot_T_n dot_i_n dot_v_n C_t_n R_t_n R_n L_n C_n T_e u_n T_d

X_n = [T_n; i_n; v_n];
X_n_dot = [dot_T_n; dot_i_n; dot_v_n];

P_heating = v_n^2/R_n;
P_loss = (T_n - T_e)/R_t_n;
dot_T_n = (P_heating-P_loss)/C_t_n;
dot_i_n = (u_n-v_n)/L_n;
dot_v_n = (i_n-v_n/R_n)/C_n;
X_n_dot = [dot_T_n; dot_i_n; dot_v_n];

v_n_equation = sqrt(R_n/R_t_n*(T_d-T_e));
i_n_equation = v_n_equation/R_n;
u_n_equation = v_n_equation;
parameters_n = [C_t_n; R_t_n; R_n; L_n; C_n; T_e; T_d];
parameters_values_n = [15.0; 3.25; 15.0; 0.004; 1.0e-5; 20.0; 250.0];
X_n_dot_sim = subs(X_n_dot, parameters_n, parameters_values_n);
u_n_sim = double(subs(u_n_equation, parameters_n, parameters_values_n));

X_n_0 = [T_e; 0; 0];
X_n_0_sim = double(subs(X_n_0, parameters_n, parameters_values_n));
tau_n = double(subs(R_t_n*C_t_n, parameters_n, parameters_values_n));
t_preheat = 7*tau_n;
t_print = 2*pi;
t_end = t_preheat+t_print;
dt_n = 1.0e-5;
t_n_sim = 0:dt_n:t_end;
n_steps = length(t_n_sim);
X_n_sim = zeros(3,n_steps);
X_n_sim(:,1) = X_n_0_sim;
X_n_dot_function = matlabFunction(subs(X_n_dot_sim,u_n,u_n_sim), 'Vars',{T_n,i_n,v_n});

for j = 1:n_steps-1
    X_n_dot_sim = X_n_dot_function(X_n_sim(1,j), X_n_sim(2,j), X_n_sim(3,j));
    X_n_sim(:,j+1) = X_n_sim(:,j) + dt_n*X_n_dot_sim;
end

figure;
plot(t_n_sim, X_n_sim(1,:), 'LineWidth', 1.5);
hold on;
yline(250,'k--','Desired temperature');
xline(t_preheat,'r--','Printing starts');
grid on;
xlabel('Time [s]');
ylabel('Nozzle temperature T_n (^{\circ}C)');
title('Q1.5: Nozzle temperature');
legend('Simulated temperature','Desired temperature', 'Printing starts','Location','best');

figure;
tiledlayout(2,1);

nexttile;
plot(t_n_sim, X_n_sim(2,:), 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Inductor current I_n (A)');
title('Q1.5: Nozzle current');

nexttile;
plot(t_n_sim, X_n_sim(3,:), 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Capacitor voltage v_n (V)');
title('Q1.5: Nozzle capacitor voltage');
%% Q1.6 Filament printing
syms s t_f rho_f c_p_f k_f w_f h_air h_bed T_air T_bed T_left T_middle T_right ds_f T_f(s,t_f) T_eq

D_f = k_f/(rho_f*c_p_f);
gamma_air = 3*h_air/(rho_f*c_p_f*w_f);
gamma_bed = h_bed/(rho_f*c_p_f*w_f);

PDE_fil = diff(T_f, t_f) == D_f*diff(T_f, s, 2) - gamma_air*(T_f - T_air) - gamma_bed*(T_f - T_bed);

T_ss = (T_left - 2*T_middle + T_right)/ds_f^2;
T_dot_f = D_f*T_ss - gamma_air*(T_middle - T_air) - gamma_bed*(T_middle - T_bed);
T_equilibrium = solve(subs(T_dot_f, [T_left, T_middle, T_right], [T_eq, T_eq, T_eq]) == 0, T_eq);

parameters_fil = [rho_f; c_p_f; k_f; w_f; h_air; h_bed; T_air; T_bed];
parameters_values_fil = [ 1250; 1800; 0.20; 0.0005; 12; 80; 20; 60];

D_f_sim = double(subs(D_f, parameters_fil, parameters_values_fil));
gamma_air_sim = double(subs(gamma_air, parameters_fil, parameters_values_fil));
gamma_bed_sim = double(subs(gamma_bed, parameters_fil, parameters_values_fil));
T_air_sim = double(subs(T_air, parameters_fil, parameters_values_fil));
T_bed_sim = double(subs(T_bed, parameters_fil, parameters_values_fil));
T_equilibrium_sim = double(subs(T_equilibrium, parameters_fil, parameters_values_fil));

a = 0.050;
b = 0.035;
omega = 1.0;
t_print = 2*pi/omega;

phi = linspace(0,2*pi,10001);
x_path = a*cos(phi);
y_path = b*sin(phi);
s_path = [0,cumsum(sqrt(diff(x_path).^2 + diff(y_path).^2))];
L_f = s_path(end);
N_f = 2410;
ds_f_sim = L_f/N_f;
s_f = (0:N_f-1)*ds_f_sim;

phi_f = interp1(s_path,phi,s_f);
t_deposit = phi_f/omega;
T_deposit = interp1(t_n_sim,X_n_sim(1,:),t_preheat);

T_dot_f_sim = subs(T_dot_f, [parameters_fil; ds_f], [parameters_values_fil; ds_f_sim]);
T_dot_f_function = matlabFunction(T_dot_f_sim, 'Vars',{T_left, T_middle, T_right});
dt_f = min(0.002, 0.4*ds_f_sim^2/D_f_sim);
tolerance_f = 1.0;
max_cooling_time = 300;
T_f_sim = nan(N_f,1);

t_now = 0;
t_history = [];
T_history = [];
while t_now < t_print + max_cooling_time
    N_active = sum(t_deposit <= t_now);
    for j = 1:N_active
        if isnan(T_f_sim(j))
            T_f_sim(j) = T_deposit;
        end
    end
    T_dot_sim = zeros(N_f,1);
    for j = 1:N_active
        if t_now >= t_print
            j_left = mod(j-2,N_f) + 1;
            j_right = mod(j,N_f) + 1;
        else
            j_left = max(j-1,1);
            j_right = min(j + 1,N_active);
        end
        T_dot_sim(j) = T_dot_f_function( T_f_sim(j_left), T_f_sim(j), T_f_sim(j_right));
    end
    T_f_sim(1:N_active) = T_f_sim(1:N_active) + dt_f*T_dot_sim(1:N_active);
    t_now = t_now + dt_f;
    if mod(round(t_now/dt_f),50) == 0
        t_history(end + 1,1) = t_now;
        T_history(end + 1,:) = T_f_sim.';
    end
    if t_now >= t_print && N_active == N_f
        T_difference = max(T_f_sim) - min(T_f_sim);

        if T_difference < tolerance_f
            break;
        end
    end
end

figure;
imagesc(s_f, t_history, T_history);
set(gca,'YDir','normal');
colorbar;
xlabel('Position along filament (m)');
ylabel('Time (s)');
title('Q1.6: Temperature distribution during printing and cooling');

figure;
plot(s_f,T_f_sim,'LineWidth',1.5);
hold on;
yline(T_equilibrium_sim,'k--','Equilibrium temperature');
grid on;
xlabel('Position along filament (m)');
ylabel('Temperature (^{\circ}C)');
title('Q1.6: Final filament temperature distribution');
legend('Filament temperature','Equilibrium temperature', 'Location','best');
%% Q3.1 Hybrid system filament replacement
r_m_sim = parameters_values(1);
R_m_sim = parameters_values(8);
K_m_sim = parameters_values(7);
L_m_sim = parameters_values(9);
C_t_n_sim = parameters_values_n(1);
R_t_n_sim = parameters_values_n(2);
R_n_sim = parameters_values_n(3);
L_n_sim = parameters_values_n(4);
C_n_sim = parameters_values_n(5);
T_e_sim = parameters_values_n(6);
T_d_sim = parameters_values_n(7);
H_sim = double(subs(H, r_m, r_m_sim));
T_low = 210;
T_high = 240;
position_tolerance = 0.0005;
velocity_tolerance = 0.003;
wipe_duration = 1.0;
load_duration = 1.0;
dt_h = 0.001;
dt_n_h = 1.0e-5;
n_nozzle = round(dt_h/dt_n_h);
t_max_h = 500;
t_move = 3.0;
t_acc = 1.0;
phi_start = 0;
a_h = parameters_values(10);
b_h = parameters_values(11);
omega_h = parameters_values(12);
x_c_h = parameters_values(13);
y_c_h = parameters_values(14);
t_print_h = 2*pi/omega_h;
xy_0_h = [x_c_h + a_h; y_c_h];
theta_0_h = H_sim\xy_0_h;
J_eq_h = double(subs(J_eq, parameters, parameters_values));
d_plus_h = double(subs(d_plus, parameters, parameters_values));
d_minus_h = double(subs(d_minus, parameters, parameters_values));

D_h = [d_plus_h, d_minus_h;
       d_minus_h, d_plus_h];
xy_dot_0_h = [0; b_h*omega_h];
xy_ddot_0_h = [-a_h*omega_h^2; 0];
theta_dot_0_h = H_sim\xy_dot_0_h;
theta_ddot_0_h = H_sim\xy_ddot_0_h;
i_0_h = (J_eq_h*theta_ddot_0_h + D_h*theta_dot_0_h)/K_m_sim;
theta_move_start = zeros(2, 1);
theta_move_end = theta_0_h(:);
X_h = [zeros(6,1); T_d_sim; u_n_sim/R_n_sim; u_n_sim];

q = "TO_FIRST_POINT";
t_h = 0;
t_print_progress = 0;
t_mode = 0;
colour = "blue";
blue_restocked = false;
colour_changed = false;
red_restocked = false;
theta_saved = theta_0_h(:);

t_history_h = [];
X_history_h = [];
q_history_h = strings(0, 1);
printed_xy = zeros(2, 0);
printed_colour = strings(1, 0);
travel_xy = zeros(2, 0);
while t_h < t_max_h && q ~= "DONE"
    theta_h = X_h(1:2);
    theta_dot_h = X_h(3:4);
    xy_h = H_sim*theta_h;
    xy_dot_h = H_sim*theta_dot_h;
    settled = norm(xy_dot_h) < velocity_tolerance;
    at_origin = norm(xy_h) < position_tolerance && settled;
    at_saved = norm(H_sim*(theta_h - theta_saved)) < position_tolerance && settled;
    q_old = q;
    switch q
        case "TO_FIRST_POINT"
            if t_h - t_mode >= t_move && at_saved
                phi_start = 0;
                q = "ACCELERATING";
            end
        case "PRINTING"
            if t_print_progress >= t_print_h
                q = "DONE";
            elseif ~blue_restocked && t_print_progress >= 0.30*t_print_h
                blue_restocked = true;
                theta_saved = theta_h(:);
                q = "STOP";
            elseif ~colour_changed && t_print_progress >= 0.50*t_print_h
                colour_changed = true;
                colour = "red";
                theta_saved = theta_h(:);
                q = "STOP";
            elseif ~red_restocked && t_print_progress >= 0.80*t_print_h
                red_restocked = true;
                theta_saved = theta_h(:);
                q = "STOP";
            end
        case "STOP"
            if settled && t_h-t_mode >= 0.2
                theta_move_start = theta_h(:);
                theta_move_end = zeros(2, 1);
                q = "TO_ORIGIN";
            end
        case "TO_ORIGIN"
            if t_h - t_mode >= t_move
                q = "WIPING";
            end
        case "WIPING"
            if t_h - t_mode >= wipe_duration
                q = "COOLING";
            end
        case "COOLING"
            if X_h(7) <= T_low
                q = "LOADING";
            end
        case "LOADING"
            if t_h - t_mode >= load_duration
                q = "HEATING";
            end
        case "HEATING"
            if X_h(7) >= T_high
                theta_move_start = theta_h(:);
                theta_move_end = theta_saved(:);
                q = "TO_SAVED";
            end
        case "TO_SAVED"
            if t_h - t_mode >= t_move && at_saved
                phi_start = omega_h*t_print_progress;
                q = "ACCELERATING";
            end
        case "ACCELERATING"
            if t_h - t_mode >= t_acc
                q = "PRINTING";
            end
    end
    if q ~= q_old
        t_mode = t_h;
    end
    if q == "DONE"
        break;
    end
    theta_d_h = theta_h(:);
    theta_dot_d_h = zeros(2, 1);
    theta_ddot_d_h = zeros(2, 1);
    theta_dddot_d_h = zeros(2, 1);
    u_nozzle_h = u_n_sim;
    switch q
        case "ACCELERATING"
            tau = min((t_h - t_mode)/t_acc,1);
            phi_h = phi_start + omega_h*t_acc*(tau^3 - 0.5*tau^4);
            phi_dot_h = omega_h*(3*tau^2-2*tau^3);
            phi_ddot_h = omega_h*(6*tau-6*tau^2)/t_acc;
            phi_dddot_h = omega_h*(6-12*tau)/t_acc^2;
            xy_d_h = [x_c_h + a_h*cos(phi_h); y_c_h + b_h*sin(phi_h)];
            xy_dot_d_h = [-a_h*sin(phi_h)*phi_dot_h; b_h*cos(phi_h)*phi_dot_h];
            xy_ddot_d_h = [-a_h*(cos(phi_h)*phi_dot_h^2 + sin(phi_h)*phi_ddot_h);
                 b_h*(-sin(phi_h)*phi_dot_h^2 + cos(phi_h)*phi_ddot_h)];
            xy_dddot_d_h = [a_h*sin(phi_h)*phi_dot_h^3 - 3*a_h*cos(phi_h)*phi_dot_h*phi_ddot_h - a_h*sin(phi_h)*phi_dddot_h; 
                -b_h*cos(phi_h)*phi_dot_h^3 - 3*b_h*sin(phi_h)*phi_dot_h*phi_ddot_h + b_h*cos(phi_h)*phi_dddot_h ];
            theta_d_h = H_sim\xy_d_h;
            theta_dot_d_h = H_sim\xy_dot_d_h;
            theta_ddot_d_h = H_sim\xy_ddot_d_h;
            theta_dddot_d_h = H_sim\xy_dddot_d_h;
        case "PRINTING"
            phi_h = omega_h*t_print_progress;
            xy_d_h = [x_c_h + a_h*cos(phi_h); y_c_h + b_h*sin(phi_h)];
            xy_dot_d_h = [-a_h*omega_h*sin(phi_h); b_h*omega_h*cos(phi_h)];
            xy_ddot_d_h = [-a_h*omega_h^2*cos(phi_h); -b_h*omega_h^2*sin(phi_h)];
            xy_dddot_d_h = [a_h*omega_h^3*sin(phi_h); -b_h*omega_h^3*cos(phi_h)];
            theta_d_h = H_sim\xy_d_h;
            theta_dot_d_h = H_sim\xy_dot_d_h;
            theta_ddot_d_h = H_sim\xy_ddot_d_h;
            theta_dddot_d_h = H_sim\xy_dddot_d_h;
        case {"TO_FIRST_POINT", "TO_ORIGIN", "TO_SAVED"}
            tau = min((t_h - t_mode)/t_move,1);
            f = 10*tau^3 - 15*tau^4 + 6*tau^5;
            f_dot = (30*tau^2 - 60*tau^3 + 30*tau^4)/t_move;
            f_ddot = (60*tau - 180*tau^2 + 120*tau^3)/t_move^2;
            f_dddot = (60 - 360*tau + 360*tau^2)/t_move^3;

            theta_move_start = theta_move_start(:);
            theta_move_end = theta_move_end(:);
            delta_theta = theta_move_end - theta_move_start;
            theta_d_h = theta_move_start + delta_theta*f;
            theta_dot_d_h = delta_theta*f_dot;
            theta_ddot_d_h = delta_theta*f_ddot;
            theta_dddot_d_h = delta_theta*f_dddot;
            if q == "TO_ORIGIN"
                u_nozzle_h = 0;
            end
        case {"STOP","WIPING", "COOLING", "LOADING", "HEATING"}
            theta_dot_d_h = zeros(2, 1);
            theta_ddot_d_h = zeros(2, 1);
            theta_dddot_d_h = zeros(2, 1);
            if q == "WIPING" || q == "COOLING" || q == "LOADING"
                u_nozzle_h = 0;
            end
    end
    theta_dot_d_h = theta_dot_d_h(:);
    theta_ddot_d_h = theta_ddot_d_h(:);
    theta_dddot_d_h = theta_dddot_d_h(:);
    i_d_h = (J_eq_h*theta_ddot_d_h + D_h*theta_dot_d_h)/K_m_sim;
    i_dot_d_h = (J_eq_h*theta_dddot_d_h + D_h*theta_ddot_d_h)/K_m_sim;
    U_h = R_m_sim*i_d_h + L_m_sim*i_dot_d_h + K_m_sim*theta_dot_d_h;
    U_h = U_h(:);

    X_dot_carriage = A_sim*X_h(1:6) + B_sim*U_h;
    X_h(1:6) = X_h(1:6) + dt_h*X_dot_carriage;
    for k = 1:n_nozzle
        T_n_h = X_h(7);
        i_n_h = X_h(8);
        v_n_h = X_h(9);
        T_n_dot_h = (v_n_h^2/R_n_sim - (T_n_h - T_e_sim)/R_t_n_sim)/C_t_n_sim;
        i_n_dot_h = (u_nozzle_h - v_n_h)/L_n_sim;
        v_n_dot_h = (i_n_h - v_n_h/R_n_sim)/C_n_sim;
        X_h(7:9) = X_h(7:9) + dt_n_h*[T_n_dot_h; i_n_dot_h; v_n_dot_h];
    end
    t_h = t_h + dt_h;
    if q ~= "PRINTING" && q ~= "ACCELERATING"
        travel_xy(:, end+1) = H_sim*X_h(1:2);
    else
        travel_xy(:, end+1) = [NaN;NaN];
    end
    if q == "PRINTING" || q == "ACCELERATING"
        if q == "ACCELERATING"
            t_print_progress = phi_h/omega_h;
        else
            t_print_progress = t_print_progress + dt_h;
        end
        printed_xy(:, end+1) = H_sim*X_h(1:2);
        printed_colour(end+1) = colour;
    end
    if mod(round(t_h/dt_h), 10) == 0
        t_history_h(end+1, 1) = t_h;
        X_history_h(end+1, :) = X_h.';
        q_history_h(end+1, 1) = q;
    end
end
figure;
hold on;
plot(travel_xy(1, :), travel_xy(2, :), 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
blue = printed_colour == "blue";
red = printed_colour == "red";
scatter(printed_xy(1, blue),printed_xy(2, blue), 5, 'b', 'filled');
scatter(printed_xy(1, red),printed_xy(2, red), 5, 'r', 'filled');
phi_plot = linspace(0, 2*pi, 501);
x_ellipse = x_c_h + a_h*cos(phi_plot);
y_ellipse = y_c_h + b_h*sin(phi_plot);
plot(x_ellipse, y_ellipse, 'k--', 'LineWidth', 1.5);
plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
axis equal;
grid on;
xlabel('x (m)');
ylabel('y (m)');
xlim([-0.05 0.25]);
ylim([-0.05 0.25]);
title('Q3.1: Hybrid printing');
legend('Non-printing movement', 'Blue filament', 'Red filament', 'Desired ellipse', 'Origin', 'Location', 'southeast');

figure;
tiledlayout(2,1);
nexttile;
plot(t_history_h, X_history_h(:, 7), 'LineWidth', 1.5);
hold on;
yline(T_low, 'r--', 'Low temperature');
yline(T_high, 'g--', 'High temperature');
grid on;
xlabel('Time (s)');
ylabel('Nozzle temperature (^{\circ}C)');
title('Q3.1: Nozzle temperature');
nexttile;
states = ["TO_FIRST_POINT", "PRINTING", "STOP", "TO_ORIGIN", "WIPING", "COOLING", "LOADING", "HEATING", "TO_SAVED", "ACCELERATING", "DONE"];
state_index = zeros(length(q_history_h), 1);
for j = 1:length(states)
    state_index(q_history_h == states(j)) = j;
end
stairs(t_history_h, state_index, 'LineWidth', 1.5);
yticks(1:length(states));
yticklabels(states);
grid on;
xlabel('Time (s)');
ylabel('FSM state');
title('Q3.1: FSM transitions');


%% Q3.2 Hybrid system - Square trajectory
square_length = 0.070;
square_half = square_length/2;
x_c_square = parameters_values(13);
y_c_square = parameters_values(14);
corners_square = [x_c_square  -  square_half, x_c_square + square_half, x_c_square + square_half, x_c_square - square_half, x_c_square - square_half;
    y_c_square - square_half, y_c_square - square_half, y_c_square + square_half, y_c_square + square_half, y_c_square - square_half];
theta_corners_square = H_sim\corners_square;
dt_square = 0.001;
dt_nozzle_square = 1.0e-5;
n_nozzle_square = round(dt_square/dt_nozzle_square);
t_side_square = 3.0;
t_move_square = 3.0;
t_max_square = 500;
restock_duration = 2.0;
T_low_square = 210;
T_high_square = 240;
angle_tol_square = 0.05;
speed_tol_square = 0.1;

X_square = [zeros(6, 1); T_d_sim; u_n_sim/R_n_sim; u_n_sim];
q_square = "TO_FIRST_CORNER";
t_square = 0;
t_mode_square = 0;
side_square = 1;
colour_square = "blue";
colour_changed_square = false;
theta_move_start_square = zeros(2, 1);
theta_move_end_square = theta_corners_square(:, 1);
theta_saved_square = theta_corners_square(:, 1);

t_history_square = [];
X_history_square = [];
q_history_square = strings(0, 1);
printed_xy_square = zeros(2, 0);
printed_colour_square = strings(1, 0);
travel_xy_square = zeros(2, 0);

while t_square < t_max_square && q_square ~= "DONE"
    theta_square = X_square(1:2);
    theta_dot_square = X_square(3:4);
    xy_square = H_sim*theta_square;
    xy_dot_square = H_sim*theta_dot_square;
    settled_square = norm(xy_dot_square) < velocity_tolerance;
    at_target_square = norm(H_sim*(theta_square - theta_move_end_square)) < position_tolerance && settled_square;
    q_old_square = q_square;
    switch q_square
        case "TO_FIRST_CORNER"
            if t_square-t_mode_square >= t_move_square && at_target_square
                q_square = "PRINTING";
            end
        case "PRINTING"
            if t_square-t_mode_square >= t_side_square
                if side_square == 4
                    q_square = "DONE";
                elseif side_square == 2 && ~colour_changed_square
                    colour_changed_square = true;
                    theta_saved_square = theta_square(:);
                    q_square = "STOP";
                else
                    side_square = side_square + 1;
                    t_mode_square = t_square;
                end
            end
        case "STOP"
            if settled_square && t_square - t_mode_square >= 0.2
                theta_move_start_square = theta_square(:);
                theta_move_end_square = zeros(2, 1);
                q_square = "TO_ORIGIN";
            end
        case "TO_ORIGIN"
            if t_square - t_mode_square >= t_move_square ...
                    && at_target_square
                q_square = "COOLING";
            end
        case "COOLING"
            if X_square(7) <= T_low_square
                q_square = "CHANGING_FILAMENT";
            end
        case "CHANGING_FILAMENT"
            if t_square - t_mode_square >= restock_duration
                colour_square = "red";
                q_square = "HEATING";
            end
        case "HEATING"
            if X_square(7) >= T_high_square
                theta_move_start_square = theta_square(:);
                theta_move_end_square = theta_saved_square(:);
                q_square = "TO_SAVED";
            end
        case "TO_SAVED"
            if t_square - t_mode_square >= t_move_square && at_target_square
                side_square = 3;
                q_square = "PRINTING";
            end
    end
    if q_square ~= q_old_square
        t_mode_square = t_square;
    end
    if q_square == "DONE"
        break;
    end

    theta_d_square = theta_square(:);
    theta_dot_d_square = zeros(2, 1);
    theta_ddot_d_square = zeros(2, 1);
    theta_dddot_d_square = zeros(2, 1);
    u_nozzle_square = u_n_sim;
    switch q_square
        case "PRINTING"
            tau = min((t_square-t_mode_square)/t_side_square,1);
            f = 10*tau^3 - 15*tau^4 + 6*tau^5;
            f_dot = (30*tau^2 - 60*tau^3 + 30*tau^4)/t_side_square;
            f_ddot = (60*tau - 180*tau^2 + 120*tau^3)/t_side_square^2;
            f_dddot = (60 - 360*tau + 360*tau^2)/t_side_square^3;
            theta_start_square = theta_corners_square(:, side_square);
            theta_end_square = theta_corners_square(:, side_square+1);
            delta_theta_square = theta_end_square - theta_start_square;
            theta_d_square = theta_start_square + delta_theta_square*f;
            theta_dot_d_square = delta_theta_square*f_dot;
            theta_ddot_d_square = delta_theta_square*f_ddot;
            theta_dddot_d_square = delta_theta_square*f_dddot;
        case {"TO_FIRST_CORNER", "TO_ORIGIN", "TO_SAVED"}
            tau = min((t_square - t_mode_square)/t_move_square, 1);
            f = 10*tau^3 - 15*tau^4 + 6*tau^5;
            f_dot = (30*tau^2 - 60*tau^3 + 30*tau^4)/t_move_square;
            f_ddot = (60*tau - 180*tau^2 + 120*tau^3)/t_move_square^2;
            f_dddot = (60 - 360*tau + 360*tau^2)/t_move_square^3;
            delta_theta_square = theta_move_end_square(:) - theta_move_start_square(:);
            theta_d_square = theta_move_start_square(:) + delta_theta_square*f;

            theta_dot_d_square = delta_theta_square*f_dot;
            theta_ddot_d_square = delta_theta_square*f_ddot;
            theta_dddot_d_square = delta_theta_square*f_dddot;
        case {"STOP","COOLING", "CHANGING_FILAMENT", "HEATING"}
            theta_dot_d_square = zeros(2, 1);
            theta_ddot_d_square = zeros(2, 1);
            theta_dddot_d_square = zeros(2, 1);
    end
    if q_square == "TO_ORIGIN" || q_square == "COOLING" || q_square == "CHANGING_FILAMENT"
        u_nozzle_square = 0;
    end
    theta_dot_d_square = theta_dot_d_square(:);
    theta_ddot_d_square = theta_ddot_d_square(:);
    theta_dddot_d_square = theta_dddot_d_square(:);
    i_d_square = (J_eq_h*theta_ddot_d_square + D_h*theta_dot_d_square)/K_m_sim;
    i_dot_d_square = (J_eq_h*theta_dddot_d_square + D_h*theta_ddot_d_square)/K_m_sim;
    U_square = R_m_sim*i_d_square + L_m_sim*i_dot_d_square + K_m_sim*theta_dot_d_square;
    U_square = U_square(:);

    X_dot_carriage_square = A_sim*X_square(1:6) + B_sim*U_square;
    X_square(1:6) = X_square(1:6) + dt_square*X_dot_carriage_square;
    for k = 1:n_nozzle_square
        T_n_square = X_square(7);
        i_n_square = X_square(8);
        v_n_square = X_square(9);
        T_n_dot_square = (v_n_square^2/R_n_sim - (T_n_square - T_e_sim)/R_t_n_sim)/C_t_n_sim;
        i_n_dot_square = (u_nozzle_square - v_n_square)/L_n_sim;
        v_n_dot_square = (i_n_square - v_n_square/R_n_sim)/C_n_sim;
        X_square(7:9) = X_square(7:9) + dt_nozzle_square*[T_n_dot_square; i_n_dot_square; v_n_dot_square];
    end
    t_square = t_square + dt_square;
    if q_square ~= "PRINTING"
        travel_xy_square(:, end+1) = H_sim*X_square(1:2);
    else
        travel_xy_square(:, end+1) = [NaN;NaN];
    end
    if q_square == "PRINTING" && X_square(7) >= T_high_square
        printed_xy_square(:, end+1) = H_sim*X_square(1:2);
        printed_colour_square(end+1) = colour_square;
    end
    if mod(round(t_square/dt_square), 10) == 0
        t_history_square(end+1, 1) = t_square;
        X_history_square(end+1, :) = X_square.';
        q_history_square(end+1, 1) = q_square;
    end
end

figure;
hold on;
plot(travel_xy_square(1, :), travel_xy_square(2, :), 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
blue_square = printed_colour_square == "blue";
red_square = printed_colour_square == "red";
scatter(printed_xy_square(1, blue_square), printed_xy_square(2, blue_square), 5, 'b', 'filled');
scatter(printed_xy_square(1, red_square), printed_xy_square(2, red_square), 5, 'r', 'filled');
plot(corners_square(1, :), corners_square(2, :), 'k--', 'LineWidth', 1.5);
plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
axis equal;
grid on;
xlim([-0.05 0.25]);
ylim([-0.05 0.25]);
xlabel('x (m)');
ylabel('y (m)');
title('Q3.2: Hybrid square printing');
legend('Non-printing movement', 'Blue filament', 'Red filament', 'Desired square', 'Origin', 'Location', 'southeast');

figure;
tiledlayout(2, 1);
nexttile;
plot(t_history_square, X_history_square(:, 7), 'LineWidth', 1.5);
hold on;
yline(T_low_square, 'r--', 'Low temperature');
yline(T_high_square, 'g--', 'High temperature');
grid on;
xlabel('Time (s)');
ylabel('Temperature (^{\circ}C)');
title('Q3.2: Nozzle temperature');
nexttile;
states_square = ["TO_FIRST_CORNER", "PRINTING", "STOP", "TO_ORIGIN", "COOLING", "CHANGING_FILAMENT", "HEATING","TO_SAVED","DONE"];
state_index_square = zeros(length(q_history_square), 1);
for j = 1:length(states_square)
    state_index_square(q_history_square == states_square(j)) = j;
end
stairs(t_history_square, state_index_square, 'LineWidth', 1.5);
yticks(1:length(states_square));
yticklabels(states_square);
grid on;
xlabel('Time (s)');
ylabel('FSM state');
title('Q3.2: FSM transitions');