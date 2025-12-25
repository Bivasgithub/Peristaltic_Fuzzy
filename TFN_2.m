x = 0:0.01:35;
% y1 = trimf(x,[-0.08, -0.07, -0.06]);
% y2 = trimf(x,[-0.11, -0.09, -0.08]);
y3 = trimf(x,[15.9986, 20.8775, 27.1492]);
plot(x,y3)
% title('trimf, P = [0.2657 0.2831.2329]78 0.2928]')
xlabel('x')
ylabel('Membership value');
hold on;
% % % ylim([-0.05 1.05])
% % Parameters for the triangular membership function
% a =0;  % Leftmost point
% b = 5;  % Peak
% c = 10; % Rightmost point
% 
% % % Triangular membership function definition
% % trimf = @(x, params) max(0, min((x - params(1)) / (params(2) - params(1)), ...
% %                                (params(3) - x) / (params(3) - params(2))));
% 
% % Define the grid
% x = 0:0.01:100; % Range for X
% y = 0:0.01:100; % Range for Y
% [X, Y] = meshgrid(x, y);
% 
% % Compute the membership values
% Z = trimf(X, [a, b, c]) .* trimf(Y, [a, b, c]);
% 
% % Plot the 3D surface
% % figure;
% surf(X, Y, Z); % Surface plot
% % colormap('viridis'); % Set colormap
% % colorbar; % Add color bar
% xlabel('X');
% ylabel('Y');
% zlabel('Membership Value');
% title('3D Triangular Membership Function');
