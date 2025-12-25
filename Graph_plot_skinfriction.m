function T1FN5

function Cf = c1 (z, A2, A3, A4, A6, Da, Gr, eta, W_s, Q, Rd, Br, Fr, M, Re, lambda11, h, dp_by_dz)
         
    Cf= -(A2*(0.0628319*cos(6.28319*z - 0.0628319) + 0.01).*((1.0*(54000.0*A2*Da*W_s + 600.0*A2*Da*W_s*sin(6.28319*z - 0.0628319) + 600.0*A2*Da*W_s*z))./(48600.0*A2*Da + 6.0*A2*Da*sin(6.28319*z - 0.0628319).^2 + 1080.0*A2*Da*z + 1080.0*A2*Da*sin(6.28319*z - 0.0628319) + 6.0*A2*Da*z.^2 + 12.0*A2*Da*z.*sin(6.28319*z - 0.0628319)) - W_s./(0.01*z + 0.01*sin(6.28319*z - 0.0628319) + 0.9) - ((15000.0*Da*Fr*Re*W_s^2 + 15000.0*Da*Fr*Re*W_s^2*lambda11).*(0.01*z + 0.01*sin(6.28319*z - 0.0628319) + 1.0).^3)./(48600.0*A2*Da + 6.0*A2*Da*sin(6.28319*z - 0.0628319).^2 + 1080.0*A2*Da*z + 1080.0*A2*Da*sin(6.28319*z - 0.0628319) + 6.0*A2*Da*z.^2 + 12.0*A2*Da*z.*sin(6.28319*z - 0.0628319)) - ((0.01*z + 0.01*sin(6.28319*z - 0.0628319) + 1.0).^2.*(18000.0*A2*W_s + 18000.0*A2*W_s*lambda11 + 200.0*A2*W_s*z + 200.0*A2*W_s*sin(6.28319*z - 0.0628319) + 18000.0*A4*Da*Gr*sin(eta) + 200.0*A2*W_s*lambda11*sin(6.28319*z - 0.0628319) + 200.0*A2*W_s*lambda11*z + 18000.0*A3*Da*M^2*W_s - 4000.0*Da*Fr*Re*W_s^2 + 18000.0*A3*Da*M^2*W_s*lambda11 - 4000.0*Da*Fr*Re*W_s^2*lambda11 + 200.0*A3*Da*M^2*W_s*z + 18000.0*A4*Da*Gr*lambda11*sin(eta) + 200.0*A4*Da*Gr*z*sin(eta) + 200.0*A3*Da*M^2*W_s*sin(6.28319*z - 0.0628319) + 200.0*A4*Da*Gr*sin(6.28319*z - 0.0628319)*sin(eta) + 200.0*A3*Da*M^2*W_s*lambda11*sin(6.28319*z - 0.0628319) + 200.0*A4*Da*Gr*lambda11*sin(6.28319*z - 0.0628319)*sin(eta) + 200.0*A3*Da*M^2*W_s*lambda11*z + 200.0*A4*Da*Gr*lambda11*z*sin(eta)))./(48600.0*A2*Da + 6.0*A2*Da*sin(6.28319*z - 0.0628319).^2 + 1080.0*A2*Da*z + 1080.0*A2*Da*sin(6.28319*z - 0.0628319) + 6.0*A2*Da*z.^2 + 12.0*A2*Da*z.*sin(6.28319*z - 0.0628319)) + (1.0*(0.01*z + 0.01*sin(6.28319*z - 0.0628319) + 1.0).*(2700.0*A2*W_s - 24300.0*Da*dp_by_dz - 3.0*Da*dp_by_dz.*sin(6.28319*z - 0.0628319).^2 + 2700.0*A2*W_s*lambda11 + 30.0*A2*W_s*z - 24300.0*Da*dp_by_dz*lambda11 - 540.0*Da*dp_by_dz.*z + 30.0*A2*W_s*sin(6.28319*z - 0.0628319) - 540.0*Da*dp_by_dz.*sin(6.28319*z - 0.0628319) - 3.0*Da*dp_by_dz.*z.^2 + 27000.0*A4*Da*Gr*sin(eta) + 30.0*A2*W_s*lambda11*sin(6.28319*z - 0.0628319) - 540.0*Da*dp_by_dz.*lambda11.*sin(6.28319*z - 0.0628319) - 3.0*Da*dp_by_dz.*lambda11.*z.^2 - 6.0*Da*dp_by_dz.*z.*sin(6.28319*z - 0.0628319) - 3.0*Da*dp_by_dz.*lambda11.*sin(6.28319*z - 0.0628319).^2 + 30.0*A2*W_s*lambda11*z + 2700.0*A3*Da*M^2*W_s - 540.0*Da*lambda11*dp_by_dz.*z - 300.0*Da*Fr*Re*W_s^2 + 3.0*A4*Da*Gr*sin(6.28319*z - 0.0628319).^2*sin(eta) + 2700.0*A3*Da*M^2*W_s*lambda11 - 300.0*Da*Fr*Re*W_s^2*lambda11 + 30.0*A3*Da*M^2*W_s*z + 27000.0*A4*Da*Gr*lambda11*sin(eta) + 570.0*A4*Da*Gr*z*sin(eta) + 30.0*A3*Da*M^2*W_s*sin(6.28319*z - 0.0628319) - 6.0*Da*lambda11*dp_by_dz.*z.*sin(6.28319*z - 0.0628319) + 570.0*A4*Da*Gr*sin(6.28319*z - 0.0628319)*sin(eta) + 3.0*A4*Da*Gr*z.^2*sin(eta) + 30.0*A3*Da*M^2*W_s*lambda11*sin(6.28319*z - 0.0628319) + 570.0*A4*Da*Gr*lambda11*sin(6.28319*z - 0.0628319)*sin(eta) + 3.0*A4*Da*Gr*lambda11*z.^2*sin(eta) + 6.0*A4*Da*Gr*z.*sin(6.28319*z - 0.0628319)*sin(eta) + 3.0*A4*Da*Gr*lambda11*sin(6.28319*z - 0.0628319).^2*sin(eta) + 30.0*A3*Da*M^2*W_s*lambda11*z + 570.0*A4*Da*Gr*lambda11*z*sin(eta) + 6.0*A4*Da*Gr*lambda11*z.*sin(6.28319*z - 0.0628319)*sin(eta)))./(48600.0*A2*Da + 6.0*A2*Da*sin(6.28319*z - 0.0628319).^2 + 1080.0*A2*Da*z + 1080.0*A2*Da*sin(6.28319*z - 0.0628319) + 6.0*A2*Da*z.^2 + 12.0*A2*Da*z.*sin(6.28319*z - 0.0628319)) + (6.94444e-7*(2.88684e9*Da*dp_by_dz + 1.1178e9*A2*W_s + 807600.0*A2*W_s*sin(6.28319*z - 0.0628319).^2 + 5880.0*A2*W_s*sin(6.28319*z - 0.0628319).^3 + 16.0*A2*W_s*sin(6.28319*z - 0.0628319).^4 + 1944000.0*Da*dp_by_dz.*sin(6.28319*z - 0.0628319).^2 + 13680.0*Da*dp_by_dz.*sin(6.28319*z - 0.0628319).^3 + 36.0*Da*dp_by_dz.*sin(6.28319*z - 0.0628319).^4 - 1.1664e10*A2*Da*W_s + 1.1178e9*A2*W_s*lambda11 + 4.914e7*A2*W_s*z + 2.88684e9*Da*dp_by_dz*lambda11 + 1.22472e8*Da*dp_by_dz.*z + 4.914e7*A2*W_s*sin(6.28319*z - 0.0628319) + 807600.0*A2*W_s*z.^2 + 5880.0*A2*W_s*z.^3 + 16.0*A2*W_s*z.^4 + 1.22472e8*Da*dp_by_dz.*sin(6.28319*z - 0.0628319) + 1944000.0*Da*dp_by_dz.*z.^2 + 13680.0*Da*dp_by_dz.*z.^3 + 36.0*Da*dp_by_dz.*z.^4 + 96.0*A2*W_s*z.^2.*sin(6.28319*z - 0.0628319).^2 - 1.76904e9*A4*Da*Gr*sin(eta) - 2.592e8*A2*Da*W_s*sin(6.28319*z - 0.0628319) + 216.0*Da*dp_by_dz.*z.^2.*sin(6.28319*z - 0.0628319).^2 - 1440000.0*A2*Da*W_s*z.^2 + 4.914e7*A2*W_s*lambda11*sin(6.28319*z - 0.0628319) + 807600.0*A2*W_s*lambda11*z.^2 + 5880.0*A2*W_s*lambda11*z.^3 + 16.0*A2*W_s*lambda11*z.^4 + 1615200.0*A2*W_s*z.*sin(6.28319*z - 0.0628319) + 1.22472e8*Da*lambda11*dp_by_dz.*sin(6.28319*z - 0.0628319) + 1944000.0*Da*dp_by_dz.*lambda11.*z.^2 + 13680.0*Da*dp_by_dz.*lambda11.*z.^3 + 36.0*Da*dp_by_dz.*lambda11.*z.^4 + 3888000.0*Da*dp_by_dz.*z.*sin(6.28319*z - 0.0628319) - 1440000.0*A2*Da*W_s*sin(6.28319*z - 0.0628319).^2 + 807600.0*A2*W_s*lambda11*sin(6.28319*z - 0.0628319).^2 + 5880.0*A2*W_s*lambda11*sin(6.28319*z - 0.0628319).^3 + 16.0*A2*W_s*lambda11*sin(6.28319*z - 0.0628319).^4 + 17640.0*A2*W_s*z.*sin(6.28319*z - 0.0628319).^2 + 17640.0*A2*W_s*z.^2.*sin(6.28319*z - 0.0628319) + 64.0*A2*W_s*z.*sin(6.28319*z - 0.0628319).^3 + 64.0*A2*W_s*z.^3.*sin(6.28319*z - 0.0628319) + 1944000.0*Da*dp_by_dz.*lambda11.*sin(6.28319*z - 0.0628319).^2 + 13680.0*Da*dp_by_dz.*lambda11.*sin(6.28319*z - 0.0628319).^3 + 36.0*Da*dp_by_dz.*lambda11.*sin(6.28319*z - 0.0628319).^4 + 41040.0*Da*dp_by_dz.*z.*sin(6.28319*z - 0.0628319).^2 + 41040.0*Da*dp_by_dz.*z.^2.*sin(6.28319*z - 0.0628319) + 144.0*Da*dp_by_dz.*z.*sin(6.28319*z - 0.0628319).^3 + 144.0*Da*dp_by_dz.*z.^3.*sin(6.28319*z - 0.0628319) - 2.592e8*A2*Da*W_s*z + 4.914e7*A2*W_s*lambda11*z + 1.1178e9*A3*Da*M^2*W_s + 1.22472e8*Da*dp_by_dz.*lambda11.*z + 6.1587e8*Da*Fr*Re*W_s^2 + 17640.0*A2*W_s*lambda11*z.*sin(6.28319*z - 0.0628319).^2 + 17640.0*A2*W_s*lambda11*z.^2.*sin(6.28319*z - 0.0628319) + 64.0*A2*W_s*lambda11*z.*sin(6.28319*z - 0.0628319).^3 + 64.0*A2*W_s*lambda11*z.^3.*sin(6.28319*z - 0.0628319) + 807600.0*A3*Da*M^2*W_s*sin(6.28319*z - 0.0628319).^2 + 5880.0*A3*Da*M^2*W_s*sin(6.28319*z - 0.0628319).^3 + 16.0*A3*Da*M^2*W_s*sin(6.28319*z - 0.0628319).^4 + 41040.0*Da*dp_by_dz.*lambda11.*z.*sin(6.28319*z - 0.0628319).^2 + 41040.0*Da*dp_by_dz.*lambda11.*z.^2.*sin(6.28319*z - 0.0628319) + 144.0*Da*dp_by_dz.*lambda11.*z.*sin(6.28319*z - 0.0628319).^3 + 144.0*Da*dp_by_dz.*lambda11.*z.^3.*sin(6.28319*z - 0.0628319) + 447600.0*Da*Fr*Re*W_s^2*sin(6.28319*z - 0.0628319).^2 + 3280.0*Da*Fr*Re*W_s^2*sin(6.28319*z - 0.0628319).^3 + 9.0*Da*Fr*Re*W_s^2*sin(6.28319*z - 0.0628319).^4 - 1136400.0*A4*Da*Gr*sin(6.28319*z - 0.0628319).^2*sin(eta) - 7800.0*A4*Da*Gr*sin(6.28319*z - 0.0628319).^3*sin(eta) - 20.0*A4*Da*Gr*sin(6.28319*z - 0.0628319).^4*sin(eta) + 1.1178e9*A3*Da*M^2*W_s*lambda11 + 6.1587e8*Da*Fr*Re*W_s^2*lambda11 + 4.914e7*A3*Da*M^2*W_s*z + 2.712e7*Da*Fr*Re*W_s^2*z + 96.0*A2*W_s*lambda11*z.^2.*sin(6.28319*z - 0.0628319).^2 - 1.76904e9*A4*Da*Gr*lambda11*sin(eta) + 216.0*Da*dp_by_dz.*lambda11.*z.^2.*sin(6.28319*z - 0.0628319).^2 - 7.3332e7*A4*Da*Gr*z*sin(eta) - 2880000.0*A2*Da*W_s*z.*sin(6.28319*z - 0.0628319) + 1615200.0*A2*W_s*lambda11*z.*sin(6.28319*z - 0.0628319) + 4.914e7*A3*Da*M^2*W_s*sin(6.28319*z - 0.0628319) + 3888000.0*Da*dp_by_dz.*lambda11.*z.*sin(6.28319*z - 0.0628319) + 807600.0*A3*Da*M^2*W_s*z.^2 + 5880.0*A3*Da*M^2*W_s*z.^3 + 16.0*A3*Da*M^2*W_s*z.^4 + 2.712e7*Da*Fr*Re*W_s^2*sin(6.28319*z - 0.0628319) + 447600.0*Da*Fr*Re*W_s^2*z.^2 + 3280.0*Da*Fr*Re*W_s^2*z.^3 + 9.0*Da*Fr*Re*W_s^2*z.^4 - 7.3332e7*A4*Da*Gr*sin(6.28319*z - 0.0628319)*sin(eta) - 1136400.0*A4*Da*Gr*z.^2*sin(eta) - 7800.0*A4*Da*Gr*z.^3*sin(eta) - 20.0*A4*Da*Gr*z.^4*sin(eta) - 120.0*A4*Da*Gr*z.^2.*sin(6.28319*z - 0.0628319).^2*sin(eta) + 4.914e7*A3*Da*M^2*W_s*lambda11*sin(6.28319*z - 0.0628319) + 807600.0*A3*Da*M^2*W_s*lambda11*z.^2 + 5880.0*A3*Da*M^2*W_s*lambda11*z.^3 + 16.0*A3*Da*M^2*W_s*lambda11*z.^4 + 2.712e7*Da*Fr*Re*W_s^2*lambda11*sin(6.28319*z - 0.0628319) + 1615200.0*A3*Da*M^2*W_s*z.*sin(6.28319*z - 0.0628319) + 447600.0*Da*Fr*Re*W_s^2*lambda11*z.^2 + 3280.0*Da*Fr*Re*W_s^2*lambda11*z.^3 + 9.0*Da*Fr*Re*W_s^2*lambda11*z.^4 + 895200.0*Da*Fr*Re*W_s^2*z.*sin(6.28319*z - 0.0628319) - 7.3332e7*A4*Da*Gr*lambda11*sin(6.28319*z - 0.0628319)*sin(eta) - 1136400.0*A4*Da*Gr*lambda11*z.^2*sin(eta) - 7800.0*A4*Da*Gr*lambda11*z.^3*sin(eta) - 20.0*A4*Da*Gr*lambda11*z.^4*sin(eta) - 2272800.0*A4*Da*Gr*z.*sin(6.28319*z - 0.0628319)*sin(eta) + 807600.0*A3*Da*M^2*W_s*lambda11*sin(6.28319*z - 0.0628319).^2 + 5880.0*A3*Da*M^2*W_s*lambda11*sin(6.28319*z - 0.0628319).^3 + 16.0*A3*Da*M^2*W_s*lambda11*sin(6.28319*z - 0.0628319).^4 + 447600.0*Da*Fr*Re*W_s^2*lambda11*sin(6.28319*z - 0.0628319).^2 + 3280.0*Da*Fr*Re*W_s^2*lambda11*sin(6.28319*z - 0.0628319).^3 + 9.0*Da*Fr*Re*W_s^2*lambda11*sin(6.28319*z - 0.0628319).^4 + 17640.0*A3*Da*M^2*W_s*z.*sin(6.28319*z - 0.0628319).^2 + 17640.0*A3*Da*M^2*W_s*z.^2.*sin(6.28319*z - 0.0628319) + 64.0*A3*Da*M^2*W_s*z.*sin(6.28319*z - 0.0628319).^3 + 64.0*A3*Da*M^2*W_s*z.^3.*sin(6.28319*z - 0.0628319) + 9840.0*Da*Fr*Re*W_s^2*z.*sin(6.28319*z - 0.0628319).^2 + 9840.0*Da*Fr*Re*W_s^2*z.^2.*sin(6.28319*z - 0.0628319) + 36.0*Da*Fr*Re*W_s^2*z.*sin(6.28319*z - 0.0628319).^3 + 36.0*Da*Fr*Re*W_s^2*z.^3.*sin(6.28319*z - 0.0628319) - 1136400.0*A4*Da*Gr*lambda11*sin(6.28319*z - 0.0628319).^2*sin(eta) - 7800.0*A4*Da*Gr*lambda11*sin(6.28319*z - 0.0628319).^3*sin(eta) - 20.0*A4*Da*Gr*lambda11*sin(6.28319*z - 0.0628319).^4*sin(eta) - 23400.0*A4*Da*Gr*z.*sin(6.28319*z - 0.0628319).^2*sin(eta) - 23400.0*A4*Da*Gr*z.^2.*sin(6.28319*z - 0.0628319)*sin(eta) - 80.0*A4*Da*Gr*z.*sin(6.28319*z - 0.0628319).^3*sin(eta) - 80.0*A4*Da*Gr*z.^3.*sin(6.28319*z - 0.0628319)*sin(eta) + 4.914e7*A3*Da*M^2*W_s*lambda11*z + 2.712e7*Da*Fr*Re*W_s^2*lambda11*z - 7.3332e7*A4*Da*Gr*lambda11*z*sin(eta) + 96.0*A3*Da*M^2*W_s*z.^2.*sin(6.28319*z - 0.0628319).^2 + 54.0*Da*Fr*Re*W_s^2*z.^2.*sin(6.28319*z - 0.0628319).^2 + 9840.0*Da*Fr*Re*W_s^2*lambda11*z.*sin(6.28319*z - 0.0628319).^2 + 9840.0*Da*Fr*Re*W_s^2*lambda11*z.^2.*sin(6.28319*z - 0.0628319) + 36.0*Da*Fr*Re*W_s^2*lambda11*z.*sin(6.28319*z - 0.0628319).^3 + 36.0*Da*Fr*Re*W_s^2*lambda11*z.^3.*sin(6.28319*z - 0.0628319) - 23400.0*A4*Da*Gr*lambda11*z.*sin(6.28319*z - 0.0628319).^2*sin(eta) - 23400.0*A4*Da*Gr*lambda11*z.^2.*sin(6.28319*z - 0.0628319)*sin(eta) - 80.0*A4*Da*Gr*lambda11*z.*sin(6.28319*z - 0.0628319).^3*sin(eta) - 80.0*A4*Da*Gr*lambda11*z.^3.*sin(6.28319*z - 0.0628319)*sin(eta) + 96.0*A3*Da*M^2*W_s*lambda11*z.^2.*sin(6.28319*z - 0.0628319).^2 + 54.0*Da*Fr*Re*W_s^2*lambda11*z.^2.*sin(6.28319*z - 0.0628319).^2 - 120.0*A4*Da*Gr*lambda11*z.^2.*sin(6.28319*z - 0.0628319).^2*sin(eta) + 1615200.0*A3*Da*M^2*W_s*lambda11*z.*sin(6.28319*z - 0.0628319) + 895200.0*Da*Fr*Re*W_s^2*lambda11*z.*sin(6.28319*z - 0.0628319) - 2272800.0*A4*Da*Gr*lambda11*z.*sin(6.28319*z - 0.0628319)*sin(eta) + 17640.0*A3*Da*M^2*W_s*lambda11*z.*sin(6.28319*z - 0.0628319).^2 + 17640.0*A3*Da*M^2*W_s*lambda11*z.^2.*sin(6.28319*z - 0.0628319) + 64.0*A3*Da*M^2*W_s*lambda11*z.*sin(6.28319*z - 0.0628319).^3 + 64.0*A3*Da*M^2*W_s*lambda11*z.^3.*sin(6.28319*z - 0.0628319)))./(Da*(log(0.01*z + 0.01*sin(6.28319*z - 0.0628319) + 1.0) + 2.30259).*(z + sin(6.28319*z - 0.0628319) + 90.0).*(0.01*z + 0.01*sin(6.28319*z - 0.0628319) + 1.0).*(90.0*A2 + A2*z + A2*sin(6.28319*z - 0.0628319)))))./(lambda11 + 1.0);
 
end
% Define the range for z
    z_values = linspace(0, 1, 50); % Adjust the range and resolution as needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_1=400;
k_2=385;
k_f=0.492;
alpha=1;
beta=1;
phi_1=0.1+0.05*alpha+beta*0.1-alpha*beta*0.1;
phi_2=0.1+0.05*alpha+beta*0.1-alpha*beta*0.1;
k_bf=k_f*((k_1+2*k_f-2*phi_1*(k_f-k_1))/(k_1+2*k_f+phi_1*(k_f-k_1)));
A6=((k_1+2*k_f-2*phi_1*(k_f-k_1))/(k_1+2*k_f+phi_1*(k_f-...
   k_1)))*((k_2*+2*k_bf-2*phi_2*(k_bf-k_2))/(k_2+2*k_bf+phi_2*(k_bf-...
   k_2)));
rho_1=8933;
rho_2=10500;
rho_f=1063;
bita_1=16.7;
bita_2=18.7;
bita_f=1.8;
sigma_f=0.6670;
sigma_1=59600000;
sigma_2=(6.3)*10^7;
A4=(1-phi_2)*((1-phi_1)+phi_1*((rho_1*bita_1)/(rho_f*bita_f)))+...
   phi_2*((rho_2*bita_2)/(rho_f*bita_f));
sigma_bf=sigma_f*((sigma_1+2*sigma_f-2*phi_1*(sigma_f-...
         sigma_1))/(sigma_1+2*sigma_f+phi_1*(sigma_f-sigma_1)));
A3=sigma_bf*((sigma_2*+2*sigma_bf-2*phi_2*(sigma_bf-...
          sigma_2))/(sigma_2+2*sigma_bf+phi_2*(sigma_bf-sigma_2)));
A2=1/(((1-phi_1)^2.5)*((1-phi_2)^2.5));
W_s=0.03; %velocity slip;
% epsilon=0.1;  % radius of inner tube
Q=0.5; %Heat sorce
Rd=0.1;
Br=0.5;
Da=0.01;
Fr=0.01;
Re=0.01;
M=0.1;%Hartmann number;
Gr=1.1; %Grashoff number;
lambda11=0.1;  
eta=pi/2; 
E1 = 0.01;
E2 = 0.02;
E3 = 0.0;
l=0.01;
m=0.01;
a=0.01;
% Define the height function h
    h = @(z) 1 + l * z + a * sin(2 * pi * (z - m));

    % Define r values based on the maximum height
    r_values = linspace(0.1, max(h(z_values)), 50);

    % Create a meshgrid for r and z
    [R, Z] = meshgrid(r_values, z_values);
% % pressure gradient define
dp_by_dz =  @(z) 8 *pi^3 * a * ( -(E1 + E2) * cos(2 * pi * (z-m)) + (E3/(2 * pi)) * sin(2 * pi * (z-m)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the temperature values for the grid
    Cf = arrayfun(@(r, z) c1(z, A2, A3, A4, A6, Da, Gr, eta, W_s, Q, Rd, Br, Fr, M, Re, lambda11, h(z), dp_by_dz(z)), R, Z);

    % Create the contour plot
%     figure;
    mesh(R, Z, Cf); % Filled contour plot with 50 levels
    xlabel('r');
    ylabel('z');
    title('skin friction');
    colormap jet
    shading interp; % Smooth the surface
    hold on;
end