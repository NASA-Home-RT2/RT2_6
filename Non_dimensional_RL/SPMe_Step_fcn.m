function [x_cse_p_updated, x_cen_updated, x_cep_updated, x_dcse_dDs_p_updated, V_out,SOC,dVout_depsilons_p,dw2,dR1,dR2, breakFlag] = ...
    SPMe_Step_fcn(I,x_cse_p,x_cen,x_cep,x_dcse_dDs_p,Dt,params)
%Step function for SPMe model with epsilon_s_p sensitivity.
%   Function receives the current state variables, current input, and time
%   increment, and outputs the updated state variables, output voltage,
%   SOC, and epsilon_s_p sensitivities.
%   Jack Fogelquist
%   2021-05-20
%   Rui Huang edited in 2021-07-15
%   SOC output changed to positive solid surface consentration only

%% Parameters
epsilon_s_n = params.epsilon_s_n; % [-] active material volume fraction of solid negative electrode
epsilon_s_p = params.epsilon_s_p; % [-] active material volume fraction of solid positive electrode
theta_n_0 = params.theta_n_0; % [-] negative electrode surface stoichiometry when SOC = 0 (JES paper = 0.0279)
theta_n_1 = params.theta_n_1; % [-] negative electrode surface stoichiometry when SOC = 1 (JES paper = 0.9014)
theta_p_0 = params.theta_p_0; % [-] positive electrode surface stoichiometry when SOC = 0 (JES paper = 0.9084)
theta_p_1 = params.theta_p_1; % [-] positive electrode surface stoichiometry when SOC = 1   (JES paper = 0.27)
%% JES other parameter
epsilon_e_sep = params.epsilon_e_sep; % [-] active material volume fraction (porosity) of electrolyte in separator
epsilon_e_n = params.epsilon_e_n; % [-] active material volume fraction (porosity) of electrolyte in negative electrode
epsilon_e = params.epsilon_e; % [-] uniform active material volume fraction of electrolyte in electrodes (for Pade approximation)
k_n = params.k_n; % [(A/m^2)*(m^3/mol)^(1-1/(3*alpha))] reaction rate in negative electrode
k_p = params.k_p; % [(A/m^2)*(m^3/mol)^(1-1/(3*alpha))] reaction rate in positive electrode
kappa = params.kappa; % [1/Ohms-m] electrolyte conductivity
t_plus = params.t_plus; % [-] transference number
D_e = params.D_e; % [m^2/s] electrolyte diffusion coefficient (assumed constant for all phases and x-coordinates; from Capiglia et al., 1999)
D_s_p = params.D_s_p; % [m^2/s] solid-phase diffusion coefficient for positive electrode
D_s_n = params.D_s_n; % [m^2/s] solid-phase diffusion coefficient for negative electrode

%% parameters in common
L_n = params.L_n; % [m] thickness of negative electrode
L_sep = params.L_sep; % [m] thickness of separator
L_p = params.L_p; % [m] thickness of positive electrode
R_s_n = params.R_s_n; % [m] radius of solid negative electrode particle
R_s_p = params.R_s_p; % [m] radius of solid positive electrode particle
Area_n = params.Area_n; % [m^2] negative electrode current collector area
Area_p = params.Area_p; % [m^2] positive electrode current collector area
F = params.F; % [C/mol] Faraday constant
R = params.R; % [J/mol-K] universal gas constant
T = params.T; % [K] cell temperature
alpha = params.alpha; % [-] charge transfer coefficient (for positive and negative electrodes)
R_f_n = params.R_f_n; % [Ohms*m^2] resistivity of SEI layer for negative electrode
R_f_p = params.R_f_p; % [Ohms*m^2] resistivity of SEI layer for positive electrode
R_c = params.R_c; %  [Ohms*m^2] contact resistance (resistance of current collectors and wiring)
cse_n_max = params.cse_n_max; % [mol/m^3] max concentration in negative electrode
cse_p_max = params.cse_p_max; % [mol/m^3] max concentration in positive electrode

gamma = (1-t_plus)/(F*Area_n); % [-] activity coefficient
a_s_n = params.a_s_n; % [1/m] negative electrode specific interfacial surface area
a_s_p = params.a_s_p; % [1/m] positive electrode specific interfacial surface area

% Initial Conditions
ce_0 = 1000; % [mol/m^3] initial concentration of Li-ions in electrolyte

%% Limits
V_out_max = params.V_out_max;
V_out_min = params.V_out_min;
SOC_max = params.SOC_max;
SOC_min = params.SOC_min;

%% initialize incremental
dw2 = 0;
dR1 = 0;
dR2 = 0;
%% SPMe Simulation

%------- ELECTROCHEMICAL COMPUTATION -------% 

%%% Intercalation Current Density
j_Li_n = I/(Area_n*L_n)*(D_s_n/R_s_n^2*Area_n*L_n*F*epsilon_s_n); % [A/m^3] deintercalation current density in negative electrode
j_Li_p = -I/(Area_p*L_p)*(D_s_p/R_s_p^2*Area_p*L_p*F*epsilon_s_p); % [A/m^3] intercalation current density in positive electrode

%%% Electrolyte Concentration (1st-Order Padé Approximation)
% State-Space Representation
C_cep = params.C_cep;
C_cen = params.C_cen;

%%% Electrode Surface Concentration (3rd-Order Padé Approximation)
C_cse_p = [3465 420 7]; % C vector
C_dcse_dDs_p = [800415 41580 903 0];

%% Discretized state transitioning matrix
G_cep = params.G_cep;
G_cen = params.G_cen;
H_cep = params.H_cep;
H_cen = params.H_cen;

G_cse_p = params.G_cse_p;
H_cse_p = params.H_cse_p;
G_dcse_dDs_p = params.G_dcse_dDs_p;
H_dcse_dDs_p = params.H_dcse_dDs_p;

%%
for i = 1:Dt
    % Update electrolyte (note that y must be updated before x because we only know the 'current' current value)
    y_cep = C_cep*x_cep;
    x_cep_updated = G_cep*x_cep + H_cep*gamma*I;
    ce_p = y_cep + ce_0;

    y_cen = C_cen*x_cen;
    x_cen_updated = G_cen*x_cen + H_cen*gamma*I;
    ce_n = y_cen + ce_0;

    % Update electrode (note that y must be updated before x because we only know the 'current' current value)
    cse_p = C_cse_p*x_cse_p; % [mol/m^3] positive electrode surface concentration
    x_cse_p_updated = G_cse_p*x_cse_p + H_cse_p*j_Li_p; % update state variables
    
    cse_p_avg = [3465 189 1] * x_cse_p_updated;
    cse_n_avg = (theta_n_0 + (cse_p_avg/cse_p_max-theta_p_0)/(theta_p_1-theta_p_0)*(theta_n_1-theta_n_0)) * cse_n_max; 
    
    %%% Concentration Polarization Electrolyte Potential Difference
    dphi_econ = 2*R*T*(1-t_plus)/F*(1+gamma)*log(ce_p/ce_n); % [V] lithium concentration polarization electrolyte potential difference

    %%% Exchange Current
    i0_n = F*k_n*(ce_0*(cse_n_max-cse_n_avg).*cse_n_avg).^alpha; % [A/m^2] negative electrode exchange current
    i0_p = F*k_p*(ce_0*(cse_p_max-cse_p).*cse_p).^alpha; % [A/m^2] positive electrode exchange current
    
    %%% Overpotential
    eta_n = R*T/(alpha*F).*log(j_Li_n./(2*a_s_n*i0_n)+sqrt((j_Li_n./(2*a_s_n*i0_n)).^2+1)); % [V] negative electrode overpotential
    eta_p = R*T/(alpha*F).*log(j_Li_p./(2*a_s_p*i0_p)+sqrt((j_Li_p./(2*a_s_p*i0_p)).^2+1)); % [V] positive electrode overpotential
    
    %%% Open Circuit Potential
    theta_n = cse_n_avg / cse_n_max; % [-] negative electrode surface stoichiometry

    U_n = 1.9793*exp(-39.3631*theta_n)+0.2482...
           -0.0909*tanh(29.8538*theta_n-0.1234)...
           -0.04478*tanh(14.9159*(theta_n-0.2769))...
           -0.0205*tanh(30.4444*(theta_n-0.6103));
       
    theta_p = cse_p / cse_p_max; % [-] positive electrode surface stoichiometry

    U_p = -0.809*theta_p+4.4875-0.0428*tanh(18.5138*(theta_p-0.5542))-...
        17.7326*tanh(15.789*(theta_p-0.3117))+...
        17.5842*tanh(15.9308*(theta_p-0.312));

    %%% Resistance
    R_lump = R_c/Area_n + R_f_p/(Area_p*L_p*a_s_p) + R_f_n/(Area_n*L_n*a_s_n) + (L_p+L_n)/(2*Area_n*kappa*epsilon_e^1.5) + L_sep/(Area_n*kappa*epsilon_e_sep^1.5); % [Ohms] lumped resistance

    %%% Terminal Voltage
    V_out = U_p - U_n + dphi_econ + eta_p - eta_n - I*R_lump*(D_s_p/R_s_p^2*Area_p*L_p*F*epsilon_s_p); % [V] terminal voltage

    %%% SOC
    SOC_p = (theta_p - theta_p_0)/(theta_p_1 - theta_p_0); % [-] SOC from positive electrode surface stoichiometry
    SOC = SOC_p; % [-] mean SOC


    %------- SENSITIVITY COMPUTATION -------%

    %%% Open Circuit Potential Derivatives WRT Surface Concentration
    dU_dcse_p = -0.809/cse_p_max ...
        -0.0428*(18.5138/cse_p_max)*((cosh(18.5138*(theta_p-0.5542)))^(-2))...
        -17.7326*(15.789/cse_p_max)*((cosh(15.789*(theta_p-0.3117)))^(-2))...
        +17.5842*(15.9308/cse_p_max)*((cosh(15.9308*(theta_p-0.312)))^(-2));

    %%% Overpotential Derivatives WRT Surface Concentration
    deta_dcse_p = -R*T/F*sign(j_Li_p)*(1+(j_Li_p/(2*a_s_p*i0_p))^-2)^(-1/2)*(cse_p_max-2*cse_p)/(cse_p*cse_p_max-cse_p^2); % [V-m^3/mol] derivative of overpotential wrt surface concentration for positive electrode

    %%% Reistance Derivatives WRT Active Material Volume Fraction
    dRlump_depsilons_p = -R_f_p*R_s_p/(3*Area_p*L_p*epsilon_s_p^2); % [Ohms] derivative of lumped resistance wrt positive electrode active material volume fraction

    %%% Overpotential Derivatives WRT Active Material Volume Fraction
    deta_das_das_depsilons_p = -R*T/(alpha*F*epsilon_s_p)*(1+(j_Li_p/(2*a_s_p*i0_p))^-2)^(-1/2)*sign(j_Li_p); % [V] derivative of overpotential wrt active material volume fraction for positive electrode

    %%% Electrode Surface Concentration Derivatives WRT Active Material Volume Fraction
    dcse_depsilons_p = cse_p - 37950.75;

    % Update dcse_dDs_p
    x_dcse_dDs_p_updated = G_dcse_dDs_p*x_dcse_dDs_p + H_dcse_dDs_p*I; % update state variables
    dcse_dDs_p = C_dcse_dDs_p*x_dcse_dDs_p_updated; 

    %%% Terminal Voltage Sensitivity WRT Active Material Volume Fraction
    dVout_depsilons_p = -dRlump_depsilons_p*I+deta_das_das_depsilons_p+(deta_dcse_p+dU_dcse_p)*dcse_depsilons_p; % [V] derivative of terminal voltage wrt positive electrode active material volume fraction
    %%% Terminal Voltage Sensitivity WRT Dsp
    dVout_dDs_p = (deta_dcse_p+dU_dcse_p)*dcse_dDs_p;

    %% update states
    x_cse_p = x_cse_p_updated;
    x_cen = x_cen_updated;
    x_cep = x_cep_updated;
    x_dcse_dDs_p = x_dcse_dDs_p_updated;
    
    dw2 = dw2 + (dVout_depsilons_p)*(dVout_dDs_p);
    dR1 = dR1 + (dVout_depsilons_p)^2;
    dR2 = dR2 + dVout_dDs_p^2;
end

%------- STOP CONDITION -------%  
breakFlag = 0; % initialize break flag
if V_out < V_out_min || V_out > V_out_max ||SOC < SOC_min || SOC > SOC_max % if voltage or SOC or ionic concentration limit is reached
    breakFlag = 1; % raise flag to indicate that the loop that is iterating this function should break
end

end