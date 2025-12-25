function T1FN
    % Function to calculate temperature distribution using stream function
    
    % Define the function Te1
    function strm = Te1(r, z, A2, A3, A4, A6, h, M, W_s, Q, Rd, Br, Da, Re, Fr, Gr, eta, lambda11, dp_by_dz)
        % Stream function calculation
        strm = (0.2*h^5*(32000.0*A2*W_s*log(10.0*h) - 320000.0*A2*W_s*h*log(10.0*h) + 32000.0*A2*W_s*lambda11*log(10.0*h) - 320000.0*A2*W_s*h*lambda11*log(10.0*h) + 32000.0*A3*Da*M^2*W_s*log(10.0*h) + 64000.0*Da*Fr*Re*W_s^2*log(10.0*h) + 32000.0*A4*Da*Gr*log(10.0*h)*sin(eta) - 320000.0*A3*Da*M^2*W_s*h*log(10.0*h) + 32000.0*A3*Da*M^2*W_s*lambda11*log(10.0*h) + 64000.0*Da*Fr*Re*W_s^2*lambda11*log(10.0*h) - 320000.0*A4*Da*Gr*h*log(10.0*h)*sin(eta) + 32000.0*A4*Da*Gr*lambda11*log(10.0*h)*sin(eta) - 320000.0*A3*Da*M^2*W_s*h*lambda11*log(10.0*h) - 320000.0*A4*Da*Gr*h*lambda11*log(10.0*h)*sin(eta)))/(28800.0*A2*Da*log(10.0*h) - 1326299.0*A2*Da*h + 6631455.0*A2*Da*h^2 - 576000.0*A2*Da*h*log(h) + 2880000.0*A2*Da*h^2*log(h)) - (0.2*r^5*(32000.0*A2*W_s*log(10.0*h) - 320000.0*A2*W_s*h*log(10.0*h) + 32000.0*A2*W_s*lambda11*log(10.0*h) - 320000.0*A2*W_s*h*lambda11*log(10.0*h) + 32000.0*A3*Da*M^2*W_s*log(10.0*h) + 64000.0*Da*Fr*Re*W_s^2*log(10.0*h) + 32000.0*A4*Da*Gr*log(10.0*h)*sin(eta) - 320000.0*A3*Da*M^2*W_s*h*log(10.0*h) + 32000.0*A3*Da*M^2*W_s*lambda11*log(10.0*h) + 64000.0*Da*Fr*Re*W_s^2*lambda11*log(10.0*h) - 320000.0*A4*Da*Gr*h*log(10.0*h)*sin(eta) + 32000.0*A4*Da*Gr*lambda11*log(10.0*h)*sin(eta) - 320000.0*A3*Da*M^2*W_s*h*lambda11*log(10.0*h) - 320000.0*A4*Da*Gr*h*lambda11*log(10.0*h)*sin(eta)))/(28800.0*A2*Da*log(10.0*h) - 1326299.0*A2*Da*h + 6631455.0*A2*Da*h^2 - 576000.0*A2*Da*h*log(h) + 2880000.0*A2*Da*h^2*log(h)) + (0.5*h^2*(36.0*Da*dp_by_dz + 20.0*A2*W_s - 51914.5*A2*Da*W_s - 200.0*A2*W_s*h + 20.0*A2*W_s*lambda11 - 720.0*Da*dp_by_dz*h + 36.0*Da*dp_by_dz*lambda11 + 12978.6*A2*W_s*h^2 - 187469.0*A2*W_s*h^3 + 576827.0*A2*W_s*h^4 + 40.0*A2*W_s*log(h) - 259572.0*Da*dp_by_dz*h^3 + 1297866.0*Da*dp_by_dz*h^4 + 72.0*Da*dp_by_dz*log(h) - 5191455.0*A2*Da*W_s*h^2 - 16.0*A4*Da*Gr*sin(eta) + 12978.6*A2*W_s*h^2*lambda11 - 187469.0*A2*W_s*h^3*lambda11 + 576827.0*A2*W_s*h^4*lambda11 - 400.0*A2*W_s*h*log(h) + 40.0*A2*W_s*lambda11*log(h) - 259572.0*Da*dp_by_dz*h^3*lambda11 + 1297866.0*Da*dp_by_dz*h^4*lambda11 - 1440.0*Da*dp_by_dz*h*log(h) + 72.0*Da*dp_by_dz*lambda11*log(h) + 7200.0*Da*dp_by_dz*h^2*log(10.0*h) + 1038299.0*A2*Da*W_s*h - 200.0*A2*W_s*h*lambda11 - 720.0*Da*dp_by_dz*h*lambda11 + 20.0*A3*Da*M^2*W_s + 13.0*Da*Fr*Re*W_s^2 - 200.0*A3*Da*M^2*W_s*h + 20.0*A3*Da*M^2*W_s*lambda11 + 13.0*Da*Fr*Re*W_s^2*lambda11 + 520.0*A4*Da*Gr*h*sin(eta) - 16.0*A4*Da*Gr*lambda11*sin(eta) - 400.0*A2*W_s*h*lambda11*log(h) + 12978.6*A3*Da*M^2*W_s*h^2 - 187469.0*A3*Da*M^2*W_s*h^3 + 576827.0*A3*Da*M^2*W_s*h^4 - 1440.0*Da*dp_by_dz*h*lambda11*log(h) + 12978.6*Da*Fr*Re*W_s^2*h^2 - 115365.0*Da*Fr*Re*W_s^2*h^3 + 324465.0*Da*Fr*Re*W_s^2*h^4 + 40.0*A3*Da*M^2*W_s*log(h) + 7200.0*Da*dp_by_dz*h^2*lambda11*log(10.0*h) + 26.0*Da*Fr*Re*W_s^2*log(h) - 3600.0*A4*Da*Gr*h^2*sin(eta) + 72103.4*A4*Da*Gr*h^3*sin(eta) - 721034.0*A4*Da*Gr*h^4*sin(eta) - 32.0*A4*Da*Gr*log(h)*sin(eta) + 12978.6*A3*Da*M^2*W_s*h^2*lambda11 - 187469.0*A3*Da*M^2*W_s*h^3*lambda11 + 576827.0*A3*Da*M^2*W_s*h^4*lambda11 - 400.0*A3*Da*M^2*W_s*h*log(h) + 12978.6*Da*Fr*Re*W_s^2*h^2*lambda11 - 115365.0*Da*Fr*Re*W_s^2*h^3*lambda11 + 324465.0*Da*Fr*Re*W_s^2*h^4*lambda11 + 40.0*A3*Da*M^2*W_s*lambda11*log(h) + 26.0*Da*Fr*Re*W_s^2*lambda11*log(h) - 3600.0*A4*Da*Gr*h^2*lambda11*sin(eta) + 72103.4*A4*Da*Gr*h^3*lambda11*sin(eta) - 721034.0*A4*Da*Gr*h^4*lambda11*sin(eta) + 1040.0*A4*Da*Gr*h*log(h)*sin(eta) - 32.0*A4*Da*Gr*lambda11*log(h)*sin(eta) - 7200.0*A4*Da*Gr*h^2*log(h)*sin(eta) - 200.0*A3*Da*M^2*W_s*h*lambda11 + 520.0*A4*Da*Gr*h*lambda11*sin(eta) - 7200.0*A4*Da*Gr*h^2*lambda11*log(h)*sin(eta) - 400.0*A3*Da*M^2*W_s*h*lambda11*log(h) + 1040.0*A4*Da*Gr*h*lambda11*log(h)*sin(eta)))/(28800.0*A2*Da*log(10.0*h) - 1326299.0*A2*Da*h + 6631455.0*A2*Da*h^2 - 576000.0*A2*Da*h*log(h) + 2880000.0*A2*Da*h^2*log(h)) - (0.5*r^2*(36.0*Da*dp_by_dz + 20.0*A2*W_s - 51914.5*A2*Da*W_s - 200.0*A2*W_s*h + 20.0*A2*W_s*lambda11 - 720.0*Da*dp_by_dz*h + 36.0*Da*dp_by_dz*lambda11 + 12978.6*A2*W_s*h^2 - 187469.0*A2*W_s*h^3 + 576827.0*A2*W_s*h^4 + 40.0*A2*W_s*log(h) - 259572.0*Da*dp_by_dz*h^3 + 1297866.0*Da*dp_by_dz*h^4 + 72.0*Da*dp_by_dz*log(h) - 5191455.0*A2*Da*W_s*h^2 - 16.0*A4*Da*Gr*sin(eta) + 12978.6*A2*W_s*h^2*lambda11 - 187469.0*A2*W_s*h^3*lambda11 + 576827.0*A2*W_s*h^4*lambda11 - 400.0*A2*W_s*h*log(h) + 40.0*A2*W_s*lambda11*log(h) - 259572.0*Da*dp_by_dz*h^3*lambda11 + 1297866.0*Da*dp_by_dz*h^4*lambda11 - 1440.0*Da*dp_by_dz*h*log(h) + 72.0*Da*dp_by_dz*lambda11*log(h) + 7200.0*Da*dp_by_dz*h^2*log(10.0*h) + 1038299.0*A2*Da*W_s*h - 200.0*A2*W_s*h*lambda11 - 720.0*Da*dp_by_dz*h*lambda11 + 20.0*A3*Da*M^2*W_s + 13.0*Da*Fr*Re*W_s^2 - 200.0*A3*Da*M^2*W_s*h + 20.0*A3*Da*M^2*W_s*lambda11 + 13.0*Da*Fr*Re*W_s^2*lambda11 + 520.0*A4*Da*Gr*h*sin(eta) - 16.0*A4*Da*Gr*lambda11*sin(eta) - 400.0*A2*W_s*h*lambda11*log(h) + 12978.6*A3*Da*M^2*W_s*h^2 - 187469.0*A3*Da*M^2*W_s*h^3 + 576827.0*A3*Da*M^2*W_s*h^4 - 1440.0*Da*dp_by_dz*h*lambda11*log(h) + 12978.6*Da*Fr*Re*W_s^2*h^2 - 115365.0*Da*Fr*Re*W_s^2*h^3 + 324465.0*Da*Fr*Re*W_s^2*h^4 + 40.0*A3*Da*M^2*W_s*log(h) + 7200.0*Da*dp_by_dz*h^2*lambda11*log(10.0*h) + 26.0*Da*Fr*Re*W_s^2*log(h) - 3600.0*A4*Da*Gr*h^2*sin(eta) + 72103.4*A4*Da*Gr*h^3*sin(eta) - 721034.0*A4*Da*Gr*h^4*sin(eta) - 32.0*A4*Da*Gr*log(h)*sin(eta) + 12978.6*A3*Da*M^2*W_s*h^2*lambda11 - 187469.0*A3*Da*M^2*W_s*h^3*lambda11 + 576827.0*A3*Da*M^2*W_s*h^4*lambda11 - 400.0*A3*Da*M^2*W_s*h*log(h) + 12978.6*Da*Fr*Re*W_s^2*h^2*lambda11 - 115365.0*Da*Fr*Re*W_s^2*h^3*lambda11 + 324465.0*Da*Fr*Re*W_s^2*h^4*lambda11 + 40.0*A3*Da*M^2*W_s*lambda11*log(h) + 26.0*Da*Fr*Re*W_s^2*lambda11*log(h) - 3600.0*A4*Da*Gr*h^2*lambda11*sin(eta) + 72103.4*A4*Da*Gr*h^3*lambda11*sin(eta) - 721034.0*A4*Da*Gr*h^4*lambda11*sin(eta) + 1040.0*A4*Da*Gr*h*log(h)*sin(eta) - 32.0*A4*Da*Gr*lambda11*log(h)*sin(eta) - 7200.0*A4*Da*Gr*h^2*log(h)*sin(eta) - 200.0*A3*Da*M^2*W_s*h*lambda11 + 520.0*A4*Da*Gr*h*lambda11*sin(eta) - 7200.0*A4*Da*Gr*h^2*lambda11*log(h)*sin(eta) - 400.0*A3*Da*M^2*W_s*h*lambda11*log(h) + 1040.0*A4*Da*Gr*h*lambda11*log(h)*sin(eta)))/(28800.0*A2*Da*log(10.0*h) - 1326299.0*A2*Da*h + 6631455.0*A2*Da*h^2 - 576000.0*A2*Da*h*log(h) + 2880000.0*A2*Da*h^2*log(h)) - (0.25*h^4*(7200.0*A2*W_s*log(10.0*h) + 7200.0*Da*dp_by_dz*log(10.0*h) + 720000.0*Da*dp_by_dz*h^2*log(10.0*h) - 72000.0*A2*W_s*h*log(10.0*h) + 7200.0*A2*W_s*lambda11*log(10.0*h) - 144000.0*Da*dp_by_dz*h*log(10.0*h) + 7200.0*Da*dp_by_dz*lambda11*log(10.0*h) - 72000.0*A2*W_s*h*lambda11*log(10.0*h) - 144000.0*Da*dp_by_dz*h*lambda11*log(10.0*h) + 7200.0*A3*Da*M^2*W_s*log(10.0*h) + 7200.0*Da*Fr*Re*W_s^2*log(10.0*h) + 720000.0*Da*dp_by_dz*h^2*lambda11*log(10.0*h) - 720000.0*A4*Da*Gr*h^2*log(10.0*h)*sin(eta) - 72000.0*A3*Da*M^2*W_s*h*log(10.0*h) + 7200.0*A3*Da*M^2*W_s*lambda11*log(10.0*h) + 7200.0*Da*Fr*Re*W_s^2*lambda11*log(10.0*h) + 72000.0*A4*Da*Gr*h*log(10.0*h)*sin(eta) - 72000.0*A3*Da*M^2*W_s*h*lambda11*log(10.0*h) + 72000.0*A4*Da*Gr*h*lambda11*log(10.0*h)*sin(eta) - 720000.0*A4*Da*Gr*h^2*lambda11*log(10.0*h)*sin(eta)))/(28800.0*A2*Da*log(10.0*h) - 1326299.0*A2*Da*h + 6631455.0*A2*Da*h^2 - 576000.0*A2*Da*h*log(h) + 2880000.0*A2*Da*h^2*log(h)) + (0.25*r^4*(7200.0*A2*W_s*log(10.0*h) + 7200.0*Da*dp_by_dz*log(10.0*h) + 720000.0*Da*dp_by_dz*h^2*log(10.0*h) - 72000.0*A2*W_s*h*log(10.0*h) + 7200.0*A2*W_s*lambda11*log(10.0*h) - 144000.0*Da*dp_by_dz*h*log(10.0*h) + 7200.0*Da*dp_by_dz*lambda11*log(10.0*h) - 72000.0*A2*W_s*h*lambda11*log(10.0*h) - 144000.0*Da*dp_by_dz*h*lambda11*log(10.0*h) + 7200.0*A3*Da*M^2*W_s*log(10.0*h) + 7200.0*Da*Fr*Re*W_s^2*log(10.0*h) + 720000.0*Da*dp_by_dz*h^2*lambda11*log(10.0*h) - 720000.0*A4*Da*Gr*h^2*log(10.0*h)*sin(eta) - 72000.0*A3*Da*M^2*W_s*h*log(10.0*h) + 7200.0*A3*Da*M^2*W_s*lambda11*log(10.0*h) + 7200.0*Da*Fr*Re*W_s^2*lambda11*log(10.0*h) + 72000.0*A4*Da*Gr*h*log(10.0*h)*sin(eta) - 72000.0*A3*Da*M^2*W_s*h*lambda11*log(10.0*h) + 72000.0*A4*Da*Gr*h*lambda11*log(10.0*h)*sin(eta) - 720000.0*A4*Da*Gr*h^2*lambda11*log(10.0*h)*sin(eta)))/(28800.0*A2*Da*log(10.0*h) - 1326299.0*A2*Da*h + 6631455.0*A2*Da*h^2 - 576000.0*A2*Da*h*log(h) + 2880000.0*A2*Da*h^2*log(h)) - (0.166667*h^6*(180000.0*Da*Fr*Re*W_s^2*log(10.0*h) + 180000.0*Da*Fr*Re*W_s^2*lambda11*log(10.0*h)))/(28800.0*A2*Da*log(10.0*h) - 1326299.0*A2*Da*h + 6631455.0*A2*Da*h^2 - 576000.0*A2*Da*h*log(h) + 2880000.0*A2*Da*h^2*log(h)) + (0.166667*r^6*(180000.0*Da*Fr*Re*W_s^2*log(10.0*h) + 180000.0*Da*Fr*Re*W_s^2*lambda11*log(10.0*h)))/(28800.0*A2*Da*log(10.0*h) - 1326299.0*A2*Da*h + 6631455.0*A2*Da*h^2 - 576000.0*A2*Da*h*log(h) + 2880000.0*A2*Da*h^2*log(h)) - (0.5*h^2*log(h)*(36.0*Da*dp_by_dz + 20.0*A2*W_s + 14400.0*A2*Da*W_s + 20.0*A2*W_s*lambda11 - 360.0*Da*dp_by_dz*h + 36.0*Da*dp_by_dz*lambda11 - 3600.0*A2*W_s*h^2 + 16000.0*A2*W_s*h^3 - 3600.0*Da*dp_by_dz*h^2 + 36000.0*Da*dp_by_dz*h^3 - 16.0*A4*Da*Gr*sin(eta) - 3600.0*A2*W_s*h^2*lambda11 + 16000.0*A2*W_s*h^3*lambda11 - 3600.0*Da*dp_by_dz*h^2*lambda11 + 36000.0*Da*dp_by_dz*h^3*lambda11 - 144000.0*A2*Da*W_s*h - 360.0*Da*dp_by_dz*h*lambda11 + 20.0*A3*Da*M^2*W_s + 13.0*Da*Fr*Re*W_s^2 + 20.0*A3*Da*M^2*W_s*lambda11 + 130.0*Da*Fr*Re*W_s^2*h + 13.0*Da*Fr*Re*W_s^2*lambda11 + 360.0*A4*Da*Gr*h*sin(eta) - 16.0*A4*Da*Gr*lambda11*sin(eta) - 3600.0*A3*Da*M^2*W_s*h^2 + 16000.0*A3*Da*M^2*W_s*h^3 - 2300.0*Da*Fr*Re*W_s^2*h^2 + 9000.0*Da*Fr*Re*W_s^2*h^3 - 20000.0*A4*Da*Gr*h^3*sin(eta) - 3600.0*A3*Da*M^2*W_s*h^2*lambda11 + 16000.0*A3*Da*M^2*W_s*h^3*lambda11 - 2300.0*Da*Fr*Re*W_s^2*h^2*lambda11 + 9000.0*Da*Fr*Re*W_s^2*h^3*lambda11 - 20000.0*A4*Da*Gr*h^3*lambda11*sin(eta) + 130.0*Da*Fr*Re*W_s^2*h*lambda11 + 360.0*A4*Da*Gr*h*lambda11*sin(eta)))/(33157.2*A2*Da - 331572.0*A2*Da*h + 14400.0*A2*Da*log(h) - 144000.0*A2*Da*h*log(h)) + (0.5*r^2*log(r)*(36.0*Da*dp_by_dz + 20.0*A2*W_s + 14400.0*A2*Da*W_s + 20.0*A2*W_s*lambda11 - 360.0*Da*dp_by_dz*h + 36.0*Da*dp_by_dz*lambda11 - 3600.0*A2*W_s*h^2 + 16000.0*A2*W_s*h^3 - 3600.0*Da*dp_by_dz*h^2 + 36000.0*Da*dp_by_dz*h^3 - 16.0*A4*Da*Gr*sin(eta) - 3600.0*A2*W_s*h^2*lambda11 + 16000.0*A2*W_s*h^3*lambda11 - 3600.0*Da*dp_by_dz*h^2*lambda11 + 36000.0*Da*dp_by_dz*h^3*lambda11 - 144000.0*A2*Da*W_s*h - 360.0*Da*dp_by_dz*h*lambda11 + 20.0*A3*Da*M^2*W_s + 13.0*Da*Fr*Re*W_s^2 + 20.0*A3*Da*M^2*W_s*lambda11 + 130.0*Da*Fr*Re*W_s^2*h + 13.0*Da*Fr*Re*W_s^2*lambda11 + 360.0*A4*Da*Gr*h*sin(eta) - 16.0*A4*Da*Gr*lambda11*sin(eta) - 3600.0*A3*Da*M^2*W_s*h^2 + 16000.0*A3*Da*M^2*W_s*h^3 - 2300.0*Da*Fr*Re*W_s^2*h^2 + 9000.0*Da*Fr*Re*W_s^2*h^3 - 20000.0*A4*Da*Gr*h^3*sin(eta) - 3600.0*A3*Da*M^2*W_s*h^2*lambda11 + 16000.0*A3*Da*M^2*W_s*h^3*lambda11 - 2300.0*Da*Fr*Re*W_s^2*h^2*lambda11 + 9000.0*Da*Fr*Re*W_s^2*h^3*lambda11 - 20000.0*A4*Da*Gr*h^3*lambda11*sin(eta) + 130.0*Da*Fr*Re*W_s^2*h*lambda11 + 360.0*A4*Da*Gr*h*lambda11*sin(eta)))/(33157.2*A2*Da - 331572.0*A2*Da*h + 14400.0*A2*Da*log(h) - 144000.0*A2*Da*h*log(h));
 
    end

    % Define the range for z
    z_values = linspace(0, 1, 50); % Adjust the range and resolution as needed

    % Set the fixed parameters
    k_1 = 400;
    k_2 = 385;
    k_f = 0.492;
    alpha = 1;
    beta = 1;

    phi_1 = 0.1 + 0.05 * alpha + beta * 0.1 - alpha * beta * 0.1;
    phi_2 = 0.1 + 0.05 * alpha + beta * 0.1 - alpha * beta * 0.1;

    k_bf = k_f * ((k_1 + 2 * k_f - 2 * phi_1 * (k_f - k_1)) / ...
                  (k_1 + 2 * k_f + phi_1 * (k_f - k_1)));

    A6 = ((k_1 + 2 * k_f - 2 * phi_1 * (k_f - k_1)) / ...
          (k_1 + 2 * k_f + phi_1 * (k_f - k_1))) * ...
         ((k_2 + 2 * k_bf - 2 * phi_2 * (k_bf - k_2)) / ...
          (k_2 + 2 * k_bf + phi_2 * (k_bf - k_2)));

    rho_1 = 8933;
    rho_2 = 10500;
    rho_f = 1063;
    bita_1 = 16.7;
    bita_2 = 18.7;
    bita_f = 1.8;

    sigma_f = 0.6670;
    sigma_1 = 59600000;
    sigma_2 = 6.3e7;

    A4 = (1 - phi_2) * ((1 - phi_1) + phi_1 * ((rho_1 * bita_1) / (rho_f * bita_f))) + ...
         phi_2 * ((rho_2 * bita_2) / (rho_f * bita_f));

    sigma_bf = sigma_f * ((sigma_1 + 2 * sigma_f - 2 * phi_1 * (sigma_f - sigma_1)) / ...
                          (sigma_1 + 2 * sigma_f + phi_1 * (sigma_f - sigma_1)));

    A3 = sigma_bf * ((sigma_2 + 2 * sigma_bf - 2 * phi_2 * (sigma_bf - sigma_2)) / ...
                     (sigma_2 + 2 * sigma_bf + phi_2 * (sigma_bf - sigma_2)));

    A2 = 1 / (((1 - phi_1)^2.5) * ((1 - phi_2)^2.5));

    W_s=0.02; %velocity slip;
    Q=1; %Heat sorce
    Rd=0.2;
    Br=0.5;
    Da=0.02;
    Fr=0.01;
    Re=0.95;
    M=0.5;%Hartmann number;
    Gr=0.01; %Grashoff number;
    lambda11=1.5;  
    eta=pi/2; 
    l = 0.01;
    m = 0.1;
    a = 0.01;

    % Define the height function h
    h = @(z) 1 + l * z + a * sin(2 * pi * (z - m));

    % Define r values based on the maximum height
    r_values = linspace(0.1, max(h(z_values)), 50);

    % Create a meshgrid for r and z
    [R, Z] = meshgrid(r_values, z_values);

    % Define the pressure gradient
    E1 = 0.16;
    E2 = 0.01;
    E3 = 0.02;

    dp_by_dz = @(z) 8 * pi^3 * a * (-(E1 + E2) * cos(2 * pi * (z - m)) + ...
                                   (E3 / (2 * pi)) * sin(2 * pi * (z - m)));

    % Calculate the temperature values for the grid
    strm = arrayfun(@(r, z) Te1(r, z, A2, A3, A4, A6, h(z), M, W_s, Q, Rd, Br, Da, Re, Fr, Gr, eta, lambda11, dp_by_dz(z)), R, Z);

    % Create the contour plot
    figure;
    contour(R, Z, strm, 120); % Filled contour plot with 50 levels
    xlabel('r');
    ylabel('z');
    title('Stream function');
    colorbar;
    shading interp; % Smooth the surface
end
