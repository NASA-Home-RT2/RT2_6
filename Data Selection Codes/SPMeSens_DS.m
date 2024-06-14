function [x_cse_n_updated, x_cse_p_updated, x_ce_updated, x_dcse_dDs_n_updated, x_dcse_dDs_p_updated, x_dcse_depsilons_n_updated, x_dcse_depsilons_p_updated, x_dce_dDe_updated, x_dce_depsilone_updated, x_dce_depsilonesep_updated, V_out, dVout_dDs_n, dVout_dDs_p, dVout_depsilons_n, dVout_depsilons_p, dVout_dk_n, dVout_dk_p, dVout_dDe, dVout_depsilone, dVout_depsilonesep, breakFlag] = SPMeSens_DS(I,Dsp_estimate,epsilon_s_p_apriori,x_cse_n,x_cse_p,x_ce,x_dcse_dDs_n,x_dcse_dDs_p,x_dcse_depsilons_n,x_dcse_depsilons_p,x_dce_dDe,x_dce_depsilone,x_dce_depsilonesep,dt,V_out_max,V_out_min,SOC_max,SOC_min)
% Step function for SPMe model with parameter sensitivity computation.
%   Jack Fogelquist
%   2021-05-20

%   Function receives the current state variables, current input, and time
%   increment, and outputs the updated state variables, output voltage, and
%   parameter sensitivities.

%   2021-06-28: Revision with '_DS' suffix includes inputs for D_s_p and
%   epsilon_s_p (implemented in data selection algorithm to compute
%   sensitivities with estimated D_s_p and a priori epsilon_s_p).

% A priori parameter values:
D_s_p = Dsp_estimate; % [m^2/s] solid-phase diffusion coefficient for positive electrode
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
epsilon_e = mean([epsilon_e_n,epsilon_e_p]); % [-] uniform active material volume fraction of electrolyte in electrodes (for Pade approximation)

% Limits
cse_n_max = 3.6e3 * 372 * 1800 / F; % [mol/m^3] max concentration in negative electrode
cse_p_max = 3.6e3 * 274 * 5010 / F; % [mol/m^3] max concentration in positive electrode
theta_n_0 = 0.0069; % [-] negative electrode surface stoichiometry when SOC = 0
theta_n_1 = 0.6760; % [-] negative electrode surface stoichiometry when SOC = 1
theta_p_0 = 0.8228; % [-] positive electrode surface stoichiometry when SOC = 0
theta_p_1 = 0.4420; % [-] positive electrode surface stoichiometry when SOC = 1

% Initial Conditions
ce_0 = 1000; % [mol/m^3] initial concentration of Li-ions in electrolyte


%% SPMe Simulation

%------- ELECTROCHEMICAL COMPUTATION -------% 

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
x_cse_n_updated = (A_cse_n*dt+eye(size(A_cse_n)))*x_cse_n + B_cse_n*dt*j_Li_n; % update state variables
cse_p = C_cse_p*x_cse_p; % [mol/m^3] positive electrode surface concentration
x_cse_p_updated = (A_cse_p*dt+eye(size(A_cse_p)))*x_cse_p + B_cse_p*dt*j_Li_p; % update state variables

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

%%% SOC
SOC_n = (theta_n - theta_n_0)/(theta_n_1 - theta_n_0); % [-] SOC from negative electrode surface stoichiometry
SOC_p = (theta_p - theta_p_0)/(theta_p_1 - theta_p_0); % [-] SOC from positive electrode surface stoichiometry
SOC = (SOC_n + SOC_p)/2; % [-] mean SOC

%------- SENSITIVITY COMPUTATION -------% 

%%% Open Circuit Potential Derivatives WRT Surface Concentration
dU_dcse_n = -1.5*(120.0/cse_n_max)*exp(-120.0*theta_n)  ...
 +(0.0351/(0.083*cse_n_max))*((cosh((theta_n-0.286)/0.083))^(-2)) ...
 -(0.0045/(cse_n_max*0.119))*((cosh((theta_n-0.849)/0.119))^(-2)) ...
 -(0.035/(cse_n_max*0.05))*((cosh((theta_n-0.9233)/0.05))^(-2)) ...
 -(0.0147/(cse_n_max*0.034))*((cosh((theta_n-0.5)/0.034))^(-2)) ...
 -(0.102/(cse_n_max*0.142))*((cosh((theta_n-0.194)/0.142))^(-2)) ...
 -(0.022/(cse_n_max*0.0164))*((cosh((theta_n-0.9)/0.0164))^(-2)) ...
 -(0.011/(cse_n_max*0.0226))*((cosh((theta_n-0.124)/0.0226))^(-2)) ...
 +(0.0155/(cse_n_max*0.029))*((cosh((theta_n-0.105)/0.029))^(-2)); % [V-m^3/mol] derivative of open circuit potential wrt surface concentration for negative electrode
dU_dcse_p = 0.07645*(-54.4806/cse_p_max)*((cosh(30.834-54.4806*theta_p))^(-2)) ...
 +2.1581*(-50.294/cse_p_max)*((cosh(52.294-50.294*theta_p))^(-2)) ...
 +0.14169*(19.8543/cse_p_max)*((cosh(11.0923-19.8543*theta_p))^(-2)) ...
 -0.2051*(5.4888/cse_p_max)*((cosh(1.4684-5.4888*theta_p))^(-2)) ...
 -0.2531/0.1316/cse_p_max*((cosh((-theta_p+0.56478)/0.1316))^(-2)) ...
 -0.02167/0.006/cse_p_max*((cosh((theta_p-0.525)/0.006))^(-2)); % [V-m^3/mol] derivative of open circuit potential wrt surface concentration for positive electrode

%%% Overpotential Derivatives WRT Surface Concentration
deta_dcse_n = -R*T/F*sign(j_Li_n)*(1+(j_Li_n/(2*a_s_n*i0_n))^-2)^(-1/2)*(cse_n_max-2*cse_n)/(cse_n*cse_n_max-cse_n^2); % [V-m^3/mol] derivative of overpotential wrt surface concentration for negative electrode
deta_dcse_p = -R*T/F*sign(j_Li_p)*(1+(j_Li_p/(2*a_s_p*i0_p))^-2)^(-1/2)*(cse_p_max-2*cse_p)/(cse_p*cse_p_max-cse_p^2); % [V-m^3/mol] derivative of overpotential wrt surface concentration for positive electrode

%%% Electrode Surface Concentration Derivatives WRT Solid-Phase Diffusion Coefficient
% Negative Electrode: State-Space Representation
A_dcse_dDs_n = [0 1 0 0; 0 0 1 0; 0 0 0 1; -12006225*D_s_n^4/R_s_n^8 -1309770*D_s_n^3/R_s_n^6 -42651*D_s_n^2/R_s_n^4 -378*D_s_n/R_s_n^2]; % A matrix
B_dcse_dDs_n = [0 0 0 1]'; % B vector
C_dcse_dDs_n = 21/(epsilon_s_n*F*R_s_n^6)*[38115*D_s_n^2 1980*R_s_n^2*D_s_n 43*R_s_n^4 0]; % C vector
% Positive Electrode: State-Space Representation
A_dcse_dDs_p = [0 1 0 0; 0 0 1 0; 0 0 0 1; -12006225*D_s_p^4/R_s_p^8 -1309770*D_s_p^3/R_s_p^6 -42651*D_s_p^2/R_s_p^4 -378*D_s_p/R_s_p^2]; % A matrix
B_dcse_dDs_p = [0 0 0 1]'; % B vector
C_dcse_dDs_p = 21/(epsilon_s_p*F*R_s_p^6)*[38115*D_s_p^2 1980*R_s_p^2*D_s_p 43*R_s_p^4 0]; % C vector
% Update
dcse_dDs_n = C_dcse_dDs_n*x_dcse_dDs_n; % [mol-s/m^5] derivative of negative electrode surface concentration wrt solid-phase diffusion coefficient
x_dcse_dDs_n_updated = (A_dcse_dDs_n*dt+eye(size(A_dcse_dDs_n)))*x_dcse_dDs_n + B_dcse_dDs_n*dt*j_Li_n; % update state variables
dcse_dDs_p = C_dcse_dDs_p*x_dcse_dDs_p; % [mol-s/m^5] derivative of positive electrode surface concentration wrt solid-phase diffusion coefficient
x_dcse_dDs_p_updated = (A_dcse_dDs_p*dt+eye(size(A_dcse_dDs_p)))*x_dcse_dDs_p + B_dcse_dDs_p*dt*j_Li_p; % update state variables

%%% Terminal Voltage Sensitivity WRT Solid-Phase Diffusion Coefficient
dVout_dDs_n = (-dU_dcse_n-deta_dcse_n)*dcse_dDs_n; % [V-s/m^2] derivative of terminal voltage wrt negative electrode solid-phase diffusion coefficient
dVout_dDs_p = (dU_dcse_p+deta_dcse_p)*dcse_dDs_p; % [V-s/m^2] derivative of terminal voltage wrt positive electrode solid-phase diffusion coefficient

%%% Reistance Derivatives WRT Active Material Volume Fraction
dRlump_depsilons_n = -R_f_n*R_s_n/(3*Area_n*L_n*epsilon_s_n^2); % [Ohms] derivative of lumped resistance wrt negative electrode active material volume fraction
dRlump_depsilons_p = -R_f_p*R_s_p/(3*Area_p*L_p*epsilon_s_p^2); % [Ohms] derivative of lumped resistance wrt positive electrode active material volume fraction

%%% Overpotential Derivatives WRT Active Material Volume Fraction
deta_das_das_depsilons_n = -R*T/(alpha*F*epsilon_s_n)*(1+(j_Li_n/(2*a_s_n*i0_n))^-2)^(-1/2)*sign(j_Li_n); % [V] derivative of overpotential wrt active material volume fraction for negative electrode
deta_das_das_depsilons_p = -R*T/(alpha*F*epsilon_s_p)*(1+(j_Li_p/(2*a_s_p*i0_p))^-2)^(-1/2)*sign(j_Li_p); % [V] derivative of overpotential wrt active material volume fraction for positive electrode

%%% Electrode Surface Concentration Derivatives WRT Active Material Volume Fraction
% Negative Electrode: State-Space Representation
A_dcse_depsilons_n = A_cse_n; % A matrix
B_dcse_depsilons_n = B_cse_n; % B vector
C_dcse_depsilons_n = -1/epsilon_s_n*C_cse_n; % C vector
% Positive Electrode: State-Space Representation
A_dcse_depsilons_p = A_cse_p; % A matrix
B_dcse_depsilons_p = B_cse_p; % B vector
C_dcse_depsilons_p = -1/epsilon_s_p*C_cse_p; % C vector
% Update
dcse_depsilons_n = C_dcse_depsilons_n*x_dcse_depsilons_n; % [mol/m^3] derivative of negative electrode surface concentration wrt active material volume fraction 
x_dcse_depsilons_n_updated = (A_dcse_depsilons_n*dt+eye(size(A_dcse_depsilons_n)))*x_dcse_depsilons_n + B_dcse_depsilons_n*dt*j_Li_n; % update state variables
dcse_depsilons_p = C_dcse_depsilons_p*x_dcse_depsilons_p; % [mol/m^3] derivative of positive electrode surface concentration wrt active material volume fraction 
x_dcse_depsilons_p_updated = (A_dcse_depsilons_p*dt+eye(size(A_dcse_depsilons_p)))*x_dcse_depsilons_p + B_dcse_depsilons_p*dt*j_Li_p; % update state variables

%%% Terminal Voltage Sensitivity WRT Active Material Volume Fraction
dVout_depsilons_n = -dRlump_depsilons_n*I-deta_das_das_depsilons_n+(-deta_dcse_n-dU_dcse_n)*dcse_depsilons_n; % [V] derivative of terminal voltage wrt negative electrode active material volume fraction
dVout_depsilons_p = -dRlump_depsilons_p*I+deta_das_das_depsilons_p+(deta_dcse_p+dU_dcse_p)*dcse_depsilons_p; % [V] derivative of terminal voltage wrt positive electrode active material volume fraction

%%% Terminal Voltage Sensitivity WRT Reaction Rate Constant
dVout_dk_n = -(-R*T/(alpha*F*k_n)*(1+(j_Li_n/(2*a_s_n*i0_n))^-2)^(-1/2)*sign(j_Li_n)); % [V/((A/m^2)*(m^3/mol)^(1-1/(3*alpha)))] derivative of terminal voltage wrt negative electrode reaction rate constant
dVout_dk_p = -R*T/(alpha*F*k_p)*(1+(j_Li_p/(2*a_s_p*i0_p))^-2)^(-1/2)*sign(j_Li_p); % [V/((A/m^2)*(m^3/mol)^(1-1/(3*alpha)))] derivative of terminal voltage wrt positive electrode reaction rate constant

%%% Electrolyte Current Collector Concentration Derivatives WRT Electrolyte Diffusion Coefficient
% Common Terms for Convenience
a_dce_dDe = (4*D_e*epsilon_e^(1/2)*epsilon_e_sep^3 + D_e*epsilon_e^2*epsilon_e_sep^(3/2));
b_dce_dDe = (epsilon_e^2*epsilon_e_sep + 24*epsilon_e^3 + 320*epsilon_e_sep^3 + 160*epsilon_e^(3/2)*epsilon_e_sep^(3/2));
% State-Space Representation
A_dce_dDe = [0 1; -3.6864e20*a_dce_dDe^2/b_dce_dDe^2 -3.84e10*a_dce_dDe/b_dce_dDe]; % A matrix
B_dce_dDe = [0 1]'; % B vector
C_dce_dDe_n = [-4.608e15*(a_dce_dDe^2/b_dce_dDe^2)*(epsilon_e^(3/2)+4*epsilon_e_sep^(3/2))/(D_e^2*epsilon_e^(3/2)*epsilon_e_sep^(3/2))*(1-t_plus)/(Area_p*F) 0]; % negative electrode: C vector
C_dce_dDe_p = -C_dce_dDe_n; % positive electrode: C vector
% Update
dce_dDe_n = C_dce_dDe_n*x_dce_dDe; % [mol-s/m^5] derivative of electrolyte current collector concentration in negative electrode wrt electrolyte diffusion coefficient
dce_dDe_p = C_dce_dDe_p*x_dce_dDe; % [mol-s/m^5] derivative of electrolyte current collector concentration in positive electrode wrt electrolyte diffusion coefficient
x_dce_dDe_updated = (A_dce_dDe*dt+eye(size(A_dce_dDe)))*x_dce_dDe + B_dce_dDe*dt*I; % update state variables

%%% Terminal Voltage Sensitivity WRT Electrolyte Diffusion Coefficient
dVout_dDe = 2*R*T*(1-t_plus)*(1+gamma)/F*(1/ce_p*dce_dDe_p-1/ce_n*dce_dDe_n); % [V-s/m^2] derivative of terminal voltage wrt electrolyte diffusion coefficient
    
%%% Electrolyte Current Collector Concentration Derivatives WRT Electrolyte Porosity
% Common Terms for Convenience
c_dce_depsilone = 3.6864e20*a_dce_dDe^2/b_dce_dDe^2*(1-t_plus)/(Area_p*F);
b0prime_dce_depsilone = (6*epsilon_e^4.5+320*epsilon_e_sep^4.5+76*epsilon_e^3*epsilon_e_sep^1.5+240*epsilon_e^1.5*epsilon_e_sep^3+3*epsilon_e^2*epsilon_e_sep^2.5)/(3.84e14*D_e^2*epsilon_e^3*epsilon_e_sep^3*(epsilon_e^1.5+4*epsilon_e_sep^1.5));
b1prime_dce_depsilone = 3/(4e4*D_e*epsilon_e^2.5);
% State-Space Representation
A_dce_depsilone = A_dce_dDe; % A matrix
B_dce_depsilone = B_dce_dDe; % B vector
C_dce_depsilone_n = [-c_dce_depsilone*b1prime_dce_depsilone -c_dce_depsilone*b0prime_dce_depsilone]; % negative electrode: C vector
C_dce_depsilone_p = -C_dce_depsilone_n; % positive electrode: C vector
% Update
dce_depsilone_n = C_dce_depsilone_n*x_dce_depsilone; % [mol/m^3] derivative of electrolyte current collector concentration in negative electrode wrt electrolyte porosity
dce_depsilone_p = C_dce_depsilone_p*x_dce_depsilone; % [mol/m^3] derivative of electrolyte current collector concentration in positive electrode wrt electrolyte porosity
x_dce_depsilone_updated = (A_dce_depsilone*dt+eye(size(A_dce_depsilone)))*x_dce_depsilone + B_dce_depsilone*dt*I; % update state variables

%%% Concentration Polarization Electrolyte Potential Difference Derivative WRT Electrolyte Porosity
dpsi_depsilone = 2*R*T*(1-t_plus)*(1+gamma)/F*(1/ce_p*dce_depsilone_p-1/ce_n*dce_depsilone_n); % [V] derivative of concentration polarization electrolyte potential difference wrt electrolyte porosity

%%% Reistance Derivative WRT Electrolyte Porosity
dRlump_depsilone = -0.75*(L_p+L_n)/(Area_n*kappa*epsilon_e^2.5); % [Ohms] derivative of lumped resistance wrt electrolyte porosity

%%% Terminal Voltage Sensitivity WRT Electrolyte Porosity
dVout_depsilone = -dRlump_depsilone*I+dpsi_depsilone; % [V] derivative of terminal voltage wrt electrolyte porosity

%%% Electrolyte Current Collector Concentration Derivatives WRT Separator Electrolyte Porosity
% Common Terms for Convenience
c_dce_depsilonesep = c_dce_depsilone;
b0prime_dce_depsilonesep = (epsilon_e^1.5-8*epsilon_e_sep^1.5-48*epsilon_e*epsilon_e_sep^0.5)/(1.536e15*D_e^2*epsilon_e_sep^3*(epsilon_e^1.5+4*epsilon_e_sep^1.5));
b1prime_dce_depsilonesep = 1.875e-5/(D_e*epsilon_e_sep^2.5);
% State-Space Representation
A_dce_depsilonesep = A_dce_dDe; % A matrix
B_dce_depsilonesep = B_dce_dDe; % B vector
C_dce_depsilonesep_n = [-c_dce_depsilonesep*b1prime_dce_depsilonesep -c_dce_depsilonesep*b0prime_dce_depsilonesep]; % negative electrode: C vector
C_dce_depsilonesep_p = -C_dce_depsilonesep_n; % positive electrode: C vector
% Update
dce_depsilonesep_n = C_dce_depsilonesep_n*x_dce_depsilonesep; % [mol/m^3] derivative of electrolyte current collector concentration in negative electrode wrt separator electrolyte porosity
dce_depsilonesep_p = C_dce_depsilonesep_p*x_dce_depsilonesep; % [mol/m^3] derivative of electrolyte current collector concentration in positive electrode wrt separator electrolyte porosity    
x_dce_depsilonesep_updated = (A_dce_depsilonesep*dt+eye(size(A_dce_depsilonesep)))*x_dce_depsilonesep + B_dce_depsilonesep*dt*I; % update state variables

%%% Concentration Polarization Electrolyte Potential Difference Derivative WRT Separator Electrolyte Porosity
dpsi_depsilonesep = 2*R*T*(1-t_plus)*(1+gamma)/F*(1/ce_p*dce_depsilonesep_p-1/ce_n*dce_depsilonesep_n); % [V] derivative of concentration polarization electrolyte potential difference wrt separator electrolyte porosity  

%%% Reistance Derivative WRT Separator Electrolyte Porosity
dRlump_depsilonesep = -1.5*L_sep/(Area_n*kappa*epsilon_e_sep^2.5); % [Ohms] derivative of lumped resistance wrt separator electrolyte porosity

%%% Terminal Voltage Sensitivity WRT Separator Electrolyte Porosity
dVout_depsilonesep = -dRlump_depsilonesep*I+dpsi_depsilonesep; % [V] derivative of terminal voltage wrt separator electrolyte porosity 


%------- STOP CONDITION -------% 

breakFlag = 0; % initialize break flag
if V_out < V_out_min || V_out > V_out_max || SOC < SOC_min || SOC > SOC_max % if voltage or SOC limit is reached
    breakFlag = 1; % raise flag to indicate that the loop that is iterating this function should break
end


end

