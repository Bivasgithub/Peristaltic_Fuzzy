syms r z A2 A3 A4 A6 C3 C4 C5 h M W_s Q Rd Br Da Re Fr Gr eta lambda11 dp_by_dz;

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
% A3=sigma_bf*((sigma_2*+2*sigma_bf-2*phi_2*(sigma_bf-...
%           sigma_2))/(sigma_2+2*sigma_bf+phi_2*(sigma_bf-sigma_2)));
% A2=1/(((1-phi_1)^2.5)*((1-phi_2)^2.5));
% W_s=0.01; %velocity slip;
% epsilon=0.1;  % radius of inner tube
% Q=0.01; %Heat sorce
% Rd=0.01;
% Br=0.01;
% Da=0.01;
% % Fr=0.01;
% Re=0.01;
% % M=0.01;%Hartmann number;
% Gr=0.01; %Grashoff number;
% % lambda11=0.01;  
% eta=pi/2; 
% % l=0.01;
% % m=0.01;
% % a=0.01;
% % % z=input('Enter z : z=-0.2:0.001:1.2 ');
% % h=1+l*z+a*sin(2*pi*(z-m));
% %Define symbolic functions
% % theta0 = ((r - h)/(epsilon - h));  % A function of r and z
% % w0= W_s * ((r - epsilon)/(h - epsilon));  % A function of r and z
a0=(A3*M^2*Br)/(A6+Rd);
a1=Q/(A6+Rd);
a2=(a0*(W_s)^2)/(h-epsilon)^2;
a3=(2*a0*((W_s)^2)*epsilon)/(h-epsilon)^2;
a4=(a0*((W_s)^2)*epsilon^2)/(h-epsilon)^2;
a5=(1+lambda11)/A2;
a6=((1+lambda11)/Da)+((A3*M^2*(1+lambda11))/A2);
a7=((Fr*Re*(1+lambda11))/A2);
a8=(A4*Gr*sin(eta)*(1+lambda11))/A2;
%%%%%%%%%%%%%%%%%%%%%%%%%pressure gradient%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E1 = 0.01;
% E2 = 0.01;
% E3 = 0.01;
% % % pressure gradient define
% dp_by_dz =  8 *pi^3 * a * ( -(E1 + E2) * cos(2 * pi * (z-m)) + (E3/(2 * pi)) * sin(2 * pi * (z-m)));
% dp_by_dz = 0.85;


w0 = ((W_s * (r - epsilon)) / (h - epsilon));

theta0 = ( (r - h) / (epsilon - h));
% % % Velocity
w11 = - r * diff(w0, r, 2) -  diff(w0, r) + a5 * (dp_by_dz)* r + a6 * w0 * r+ ...
        a7 * w0^2 * r - a8 * theta0 * r;
% % % Indefinite integration with respect to r
w12 = int(w11, r) + C3;

% Integrate theta12 with respect to r
w1 = int((w12 / (r)), r) + C4;

% Apply the boundary conditions
eq3 = subs(w1, r, epsilon) == 0; % Applying boundary condition at r = epsilon
eq4 = subs(w1, r, h) == 0;       % Applying boundary condition at r = h

% Solve the system of equations for C1 and C2
solution = solve([eq3, eq4], [C3, C4]);

% Substitute C1 and C2 back into theta
w1 = subs(w1, [C3, C4], [solution.C3, solution.C4]);

% % Sum the symbolic functions u0 and u1
w = w0 + w1;
psi = int((r * w), r) + C5;
% Apply the boundary conditions
eq5 = subs(psi, r, h) == 0; % Applying boundary condition at r = h
% Solve the system of equations for C1 and C2
solution = solve(eq5, C5);
psi = subs(psi, C5, solution);
% Display the result
% Display the final result rounded to 6 decimal places
F_final = vpa(psi, 6);
% disp('Function w:');
% disp(F_final);
disp('Function \psi:');
disp(F_final);