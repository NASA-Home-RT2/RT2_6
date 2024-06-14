function [x_cse_n_updated, x_cse_p_updated, x_ce_updated, V_out] = SPMe(I,x_cse_n,x_cse_p,x_ce,D_s_p_hat,epsilon_s_p_apriori,dt)
% Step function for SPMe model without parameter sensitivity computation.
%   Jack Fogelquist
%   2021-04-21

%   Function receives the current state variables, current input, and time
%   increment, and outputs the updated state variables and voltage.

%   2021-06-15: A priori value of epsilon_s_p is included as an input.

% A priori parameter values:
D_s_p = D_s_p_hat; % [m^2/s] solid-phase diffusion coefficient for positive electrode
epsilon_s_p = epsilon_s_p_apriori; % [-] active material volume fraction of solid positive electrode

% Model Parameters
C_bat = 25.67; % [Ah] battery capacity
L_n = 100e-6; % [m] thickness of negative electrode
L_sep = 25e-6; % [m] thickness of separator
L_p = 100e-6; % [m] thickness of positive electrode
R_s_n = 10e-6; % [m] radius of solid negative electrode particle
R_s_p = 10e-6; % [m] radius of solid positive electrode particle
Area_n = 1; % [m^2] negative electrode current collector area
Area_p = 1; % [m^2] positive electrode current collector area
epsilon_e_n = 0.3; % [-] active material volume fraction (porosity) of electrolyte in negative electrode
epsilon_e_sep = 1.0; % [-] active material volume fraction (porosity) of electrolyte in separator
epsilon_e_p = 0.3; % [-] active material volume fraction (porosity) of electrolyte in positive electrode
epsilon_s_n = 0.6; % [-] active material volume fraction of solid negative electrode
% epsilon_s_p = 0.5; % [-] active material volume fraction of solid positive electrode
F = 96485.33289; % [C/mol] Faraday constant
t_plus = 0.4; % [-] transference number
D_e = 2.7877e-10; % [m^2/s] electrolyte diffusion coefficient (assumed constant for all phases and x-coordinates; from Capiglia et al., 1999)
D_s_n = 3.9e-14; % [m^2/s] solid-phase diffusion coefficient for negative electrode
% D_s_p = 1e-13; % [m^2/s] solid-phase diffusion coefficient for positive electrode
R = 8.314472; % [J/mol-K] universal gas constant
T = 298.15; % [K] cell temperature
gamma = 1.2383; % [-] activity coefficient
k_n = 1e-5 / F; % [(A/m^2)*(m^3/mol)^(1-1/(3*alpha))] reaction rate in negative electrode
k_p = 3e-7 / F; % [(A/m^2)*(m^3/mol)^(1-1/(3*alpha))] reaction rate in positive electrode
alpha = 0.5; % [-] charge transfer coefficient (for positive and negative electrodes)
R_f_n = 1e-3; % [Ohms*m^2] resistivity of SEI layer for negative electrode
R_f_p = 0; % [Ohms*m^2] resistivity of SEI layer for positive electrode
R_c = 0; %  [Ohms*m^2] contact resistance (resistance of current collectors and wiring)
kappa = 1.1046; % [1/Ohms-m] electrolyte conductivity

% Computed Parameters
a_s_n = 3*epsilon_s_n/R_s_n; % [1/m] negative electrode specific interfacial surface area
a_s_p = 3*epsilon_s_p/R_s_p; % [1/m] positive electrode specific interfacial surface area
epsilon_e = (epsilon_e_n+epsilon_e_p)/2; % [-] uniform active material volume fraction of electrolyte in electrodes (for Pade approximation)

% Limits
cse_n_max = 3.6e3 * 372 * 1800 / F; % [mol/m^3] max concentration in negative electrode
cse_p_max = 3.6e3 * 274 * 5010 / F; % [mol/m^3] max concentration in positive electrode

% Initial Conditions
ce_0 = 1000; % [mol/m^3] initial concentration of Li-ions in electrolyte


%% SPMe Simulation
%%% Intercalation Current Density
j_Li_n = I/(Area_n*L_n); % [A/m^3] deintercalation current density in negative electrode
j_Li_p = -I/(Area_p*L_p); % [A/m^3] intercalation current density in positive electrode

%%% Electrolyte Concentration (1st-Order Padé Approximation)
% State-Space Representation
A_ce = -(19200000000*(4*D_e*epsilon_e^(1/2)*epsilon_e_sep^3 + D_e*epsilon_e^2*epsilon_e_sep^(3/2)))/(epsilon_e^2*epsilon_e_sep + 24*epsilon_e^3 + 320*epsilon_e_sep^3 + 160*epsilon_e^(3/2)*epsilon_e_sep^(3/2)); % A coefficient
B_ce = 1; % B coefficient
C_ce = (-1)*(240000*(epsilon_e^(3/2) + 4*epsilon_e_sep^(3/2))*(4*D_e*epsilon_e^(1/2)*epsilon_e_sep^3 + D_e*epsilon_e^2*epsilon_e_sep^(3/2))*(1 - t_plus))/(Area_p*D_e*F*epsilon_e^(3/2)*epsilon_e_sep^(3/2)*(epsilon_e^2*epsilon_e_sep + 24*epsilon_e^3 + 320*epsilon_e_sep^3 + 160*epsilon_e^(3/2)*epsilon_e_sep^(3/2))); % C coefficient
% Update (note that y must be updated before x because we only know the 'current' current value)
y_ce = C_ce*x_ce; % compute positive current collector concentration
x_ce_updated = (A_ce*dt+1)*x_ce + B_ce*dt*I; % update state variable
ce_p = y_ce + ce_0; % [mol/m^3] positive current collector concentration: apply initial concentration offset
ce_n = -y_ce + ce_0; % [mol/m^3] negative current collector concentration: reverse the direction of y_ce (symmetry) and apply initial concentration offset

%%% Electrode Surface Concentration (3rd-Order Padé Approximation)
% Negative Electrode: Derived Parameters
m_n = 1/(D_s_n*a_s_n*F); % [s-mol/m-C] convenience variable
a1_cse_n = 189*D_s_n/R_s_n^2;
a2_cse_n = 3465*D_s_n^2/R_s_n^4;
b0_cse_n = -21*m_n*D_s_n/R_s_n;
b1_cse_n = -1260*m_n*D_s_n^2/R_s_n^3;
b2_cse_n = -10395*m_n*D_s_n^3/R_s_n^5;
% Negative Electrode: State-Space Representation
A_cse_n = [0 1 0; 0 0 1; 0 -a2_cse_n -a1_cse_n]; % A matrix
B_cse_n = [0 0 1]'; % B vector
C_cse_n = [b2_cse_n b1_cse_n b0_cse_n]; % C vector
% Positive Electrode: Derived Parameters
m_p = 1/(D_s_p*a_s_p*F); % [s-mol/m-C] convenience variable
a1_cse_p = 189*D_s_p/R_s_p^2;
a2_cse_p = 3465*D_s_p^2/R_s_p^4;
b0_cse_p = -21*m_p*D_s_p/R_s_p;
b1_cse_p = -1260*m_p*D_s_p^2/R_s_p^3;
b2_cse_p = -10395*m_p*D_s_p^3/R_s_p^5;
% Positive Electrode: State-Space Representation
A_cse_p = [0 1 0; 0 0 1; 0 -a2_cse_p -a1_cse_p]; % A matrix
B_cse_p = [0 0 1]'; % B vector
C_cse_p = [b2_cse_p b1_cse_p b0_cse_p]; % C vector
% Update (note that y must be updated before x because we only know the 'current' current value)
cse_n = C_cse_n*x_cse_n; % [mol/m^3] negative electrode surface concentration
x_cse_n_updated = (A_cse_n*dt+[1 0 0; 0 1 0; 0 0 1])*x_cse_n + B_cse_n*dt*j_Li_n; % update state variables
cse_p = C_cse_p*x_cse_p; % [mol/m^3] positive electrode surface concentration
x_cse_p_updated = (A_cse_p*dt+[1 0 0; 0 1 0; 0 0 1])*x_cse_p + B_cse_p*dt*j_Li_p; % update state variables

%%% Concentration Polarization Electrolyte Potential Difference
dphi_econ = 2*R*T*(1-t_plus)/F*(1+gamma)*log(ce_p/ce_n); % [V] lithium concentration polarization electrolyte potential difference

%%% Exchange Current
i0_n = F*k_n*(ce_0*(cse_n_max-cse_n)*cse_n)^alpha; % [A/m^2] negative electrode exchange current
i0_p = F*k_p*(ce_0*(cse_p_max-cse_p)*cse_p)^alpha; % [A/m^2] positive electrode exchange current

%%% Overpotential
eta_n = R*T/(alpha*F)*log(j_Li_n/(2*a_s_n*i0_n)+sqrt((j_Li_n/(2*a_s_n*i0_n))^2+1)); % [V] negative electrode overpotential
eta_p = R*T/(alpha*F)*log(j_Li_p/(2*a_s_p*i0_p)+sqrt((j_Li_p/(2*a_s_p*i0_p))^2+1)); % [V] positive electrode overpotential

%%% Open Circuit Potential
theta_n = cse_n / cse_n_max; % [-] negative electrode surface stoichiometry
U_n = 0.194+1.5*exp(-120.0*theta_n) ...
 + 0.0351*tanh((theta_n-0.286)/0.083) ... 
 - 0.0045*tanh((theta_n-0.849)/0.119) ...
 - 0.035*tanh((theta_n-0.9233)/0.05) ...
 - 0.0147*tanh((theta_n-0.5)/0.034) ...
 - 0.102*tanh((theta_n-0.194)/0.142) ...
 - 0.022*tanh((theta_n-0.9)/0.0164) ...
 - 0.011*tanh((theta_n-0.124)/0.0226) ...
 + 0.0155*tanh((theta_n-0.105)/0.029); % [V] open circuit potential for negative electrode
theta_p = cse_p / cse_p_max; % [-] positive electrode surface stoichiometry
U_p = 2.16216+0.07645*tanh(30.834-54.4806*theta_p) ...
 + 2.1581*tanh(52.294-50.294*theta_p) ...
 - 0.14169*tanh(11.0923-19.8543*theta_p) ...
 + 0.2051*tanh(1.4684-5.4888*theta_p) ...
 + 0.2531*tanh((-theta_p+0.56478)/0.1316) ...
 - 0.02167*tanh((theta_p-0.525)/0.006); % [V] open circuit potential for positive electrode

%%% Resistance
R_lump = R_c/Area_n + R_f_p/(Area_p*L_p*a_s_p) + R_f_n/(Area_n*L_n*a_s_n) + (L_p+L_n)/(2*Area_n*kappa*epsilon_e^1.5) + L_sep/(Area_n*kappa*epsilon_e_sep^1.5); % [Ohms] lumped resistance

%%% Terminal Voltage
V_out = U_p - U_n + dphi_econ + eta_p - eta_n - I*R_lump; % [V] terminal voltage

end

