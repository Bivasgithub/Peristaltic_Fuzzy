% Parameters for the triangular membership function
c = -0.05;  % Leftmost point
b = -0.06;  % Peak
a = -0.07; % Rightmost point

% % Triangular membership function definition
% trimf = @(x, params) max(0, min((x - params(1)) / (params(2) - params(1)), ...
%                                (params(3) - x) / (params(3) - params(2))));

% Define the range for X
x = -0.2:0.01:0;

% Compute the membership values
y = trimf(x, [a, b, c]);

% Plot the 2D triangular membership function
figure;
plot(x, y, 'b-', 'LineWidth', 2);
xlabel('X');
ylabel('Membership Value');
title('2D Triangular Membership Function');
grid on;
