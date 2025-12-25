syms r z A6 C1 C2 h W_s Q Rd Br M;

%%%%%%%%%%%%%%%%%%%%%%%%%%parameter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%
% k_1=400;
% k_2=385;
% k_f=0.492;
% % alpha=0.25;
% % beta=0.25;
% phi_1=0.1+0.05*alpha+beta*0.1-alpha*beta*0.1;
% phi_2=0.1+0.05*alpha+beta*0.1-alpha*beta*0.1;
% % phi_1=0.03;
% % phi_2=0.03;
% k_bf=k_f*((k_1+2*k_f-2*phi_1*(k_f-k_1))/(k_1+2*k_f+phi_1*(k_f-k_1)));
% % k_hnf=k_f*(((k_1+2*k_f-2*phi_1*(k_f-k_1))/(k_1+2*k_f+phi_1*(k_f-...
% %       k_1)))*((k_2*+2*k_bf-2*phi_2*(k_bf-k_2))/(k_2+2*k_bf+...
% %       phi_2*(k_bf-k_2))));
% A6=((k_1+2*k_f-2*phi_1*(k_f-k_1))/(k_1+2*k_f+phi_1*(k_f-...
%    k_1)))*((k_2*+2*k_bf-2*phi_2*(k_bf-k_2))/(k_2+2*k_bf+phi_2*(k_bf-...
%    k_2)));
% % x1=1/(((1-phi_1)^2.5)*((1-phi_2)^2.5)*((1-phi_3)^2.5));
% rho_1=8933;
% rho_2=10500;
% % rho_3=19320;
% rho_f=1063;
% bita_1=16.7;
% bita_2=18.7;
% % bita_3=14;
% bita_f=1.8;
% sigma_f=0.6670;
% sigma_1=59600000;
% sigma_2=(6.3)*10^7;
% % sigma_3=(4.10)*10^7;
% A4=(1-phi_2)*((1-phi_1)+phi_1*((rho_1*bita_1)/(rho_f*bita_f)))+...
%    phi_2*((rho_2*bita_2)/(rho_f*bita_f));
% sigma_bf=sigma_f*((sigma_1+2*sigma_f-2*phi_1*(sigma_f-...
%          sigma_1))/(sigma_1+2*sigma_f+phi_1*(sigma_f-sigma_1)));
A3=sigma_bf*((sigma_2*+2*sigma_bf-2*phi_2*(sigma_bf-...
%           sigma_2))/(sigma_2+2*sigma_bf+phi_2*(sigma_bf-sigma_2)));
% A2=1/(((1-phi_1)^2.5)*((1-phi_2)^2.5));
% W_s=0.01; %velocity slip;
epsilon=0.1;  % radius of inner tube
% Q=0.01; %Heat sorce
% Rd=0.01;
% Br=0.01;
% M=0.01;%Hartmann number;
l=0.01;
m=0.01;
a=0.01;
% z=input('Enter z : z=-0.2:0.001:1.2 ');
h=1+l*z+a*sin(2*pi*(z-m));
% %Define symbolic functions
% % theta0 = ((r - h)/(epsilon - h));  % A function of r and z
% % w0= W_s * ((r - epsilon)/(h - epsilon));  % A function of r and z
a0=(A3*M^2*Br)/(A6+Rd);
a1=Q/(A6+Rd);
w0 = ((W_s * (r - epsilon)) / (h - epsilon));
theta0 = ( (r - h) / (epsilon - h));
% % % Velocity
theta11 = - r * diff(theta0, r, 2) -  diff(theta0, r)- a0 * w0^2 * r - a1 * r;
        
% % % Indefinite integration with respect to r
theta12 = int(theta11, r) + C1;

% Integrate theta12 with respect to r
theta1 = int((theta12 / (r)), r) + C2;

% Apply the boundary conditions
eq3 = subs(theta1, r, epsilon) == 0; % Applying boundary condition at r = epsilon
eq4 = subs(theta1, r, h) == 0;       % Applying boundary condition at r = h

% Solve the system of equations for C1 and C2
solution = solve([eq3, eq4], [C1, C2]);

% Substitute C1 and C2 back into theta
theta1 = subs(theta1, [C1, C2], [solution.C1, solution.C2]);

% % Sum the symbolic functions u0 and u1
theta = theta0 + theta1;
Nu =  A6 * diff(h, z) * (diff(theta, r));
Nu_at_h = subs(Nu, r, h);
% Display the result
% Display the final result rounded to 6 decimal places
F_final = vpa(Nu_at_h, 6);
% disp('Function w:');
% disp(F_final);
disp('Function Nu:');
disp(F_final);