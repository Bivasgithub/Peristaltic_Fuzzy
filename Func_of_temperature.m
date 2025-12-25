% Clear all previous variables
clear all;
clc;

% Define symbolic variables
syms r z alpha beta C1 C2 h;

%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_1 = 400;
k_2 = 385;
k_f = 0.492;

% Volume fractions for the phases
phi_1 = 0.1 + 0.05 * alpha + beta * 0.1 - alpha * beta * 0.1;
phi_2 = 0.1 + 0.05 * alpha + beta * 0.1 - alpha * beta * 0.1;

% Effective thermal conductivity for base fluid
k_bf = k_f * ((k_1 + 2 * k_f - 2 * phi_1 * (k_f - k_1)) / (k_1 + 2 * k_f + phi_1 * (k_f - k_1)));

% Effective thermal conductivity for the composite
A6 = ((k_1 + 2 * k_f - 2 * phi_1 * (k_f - k_1)) / (k_1 + 2 * k_f + phi_1 * (k_f - k_1))) * ...
     ((k_2 + 2 * k_bf - 2 * phi_2 * (k_bf - k_2)) / (k_2 + 2 * k_bf + phi_2 * (k_bf - k_2)));

% Parameters for density and thermal expansion
rho_1 = 8933;
rho_2 = 10500;
rho_f = 1063;
beta_1 = 16.7;
beta_2 = 18.7;
beta_f = 1.8;

% Electrical conductivity
sigma_f = 0.6670;
sigma_1 = 59600000;
sigma_2 = 6.3 * 10^7;

A4 = (1 - phi_2) * ((1 - phi_1) + phi_1 * ((rho_1 * beta_1) / (rho_f * beta_f))) + ...
     phi_2 * ((rho_2 * beta_2) / (rho_f * beta_f));

sigma_bf = sigma_f * ((sigma_1 + 2 * sigma_f - 2 * phi_1 * (sigma_f - ...
           sigma_1)) / (sigma_1 + 2 * sigma_f + phi_1 * (sigma_f - sigma_1)));
A3 = sigma_bf * ((sigma_2 + 2 * sigma_bf - 2 * phi_2 * (sigma_bf -...
     sigma_2)) / (sigma_2 + 2 * sigma_bf + phi_2 * (sigma_bf - sigma_2)));

A2 = 1 / (((1 - phi_1)^2.5) * ((1 - phi_2)^2.5));

% Constants
epsilon = 0.1;  % Radius of inner tube
Q = 3;        % Heat source
Rd = 0.5;
l = 0.01;
m = 0.01;
a = 0.01;
h = 1 + l * z + a * sin(2 * pi * (z - m)); % Definition of h
% Q = 0.1;
% Rd = 0.5;
theta0 = ((r - h) / (epsilon - h));
% Calculation
a1 = Q/(A6 + Rd);
% theta11 = -r * Q/(A6 + Rd);
theta11 = - r * diff(theta0, r, 2) -  diff(theta0, r) - a1 * r;
% Indefinite integration with respect to r
theta12 = int(theta11, r) + C1; 

% Integrate theta12 with respect to r
theta = int((theta12 / r), r) + C2; 

% Apply the boundary conditions
eq1 = subs(theta, r, epsilon) == 1; % Applying boundary condition at r = epsilon
eq2 = subs(theta, r, h) == 0;       % Applying boundary condition at r = h

% Solve the system of equations for C1 and C2
solution = solve([eq1, eq2], [C1, C2]);

% Substitute C1 and C2 back into theta
theta = subs(theta, [C1, C2], [solution.C1, solution.C2]);

% Display the final result rounded to 6 decimal places
theta_final = vpa(theta, 2);
% disp('The integrated function with the boundary conditions is: ');
disp(theta);
