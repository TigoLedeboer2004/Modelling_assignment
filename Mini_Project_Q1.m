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
Z = [xy_matrix(1); xy_dot_matrix(1); i_x; xy_matrix(2); xy_dot_matrix(2); i_y]
T = [r_m/2, r_m/2, 0, 0, 0, 0;
    0, 0, r_m/2, r_m/2, 0, 0;
    0, 0, 0, 0, 1, 1;
    r_m/2, r_m/2, 0, 0, 0, 0;
    0, 0, -r_m/2, -r_m/2, 0, 0;
    0, 0, 0, 0, 1, -1];

simplify(Z - T*X)