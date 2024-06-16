function [x_cse_n_updated, x_cse_p_updated, x_cen_updated, x_cep_updated, x_dcse_depsilons_p_updated, V_out, SOC, dVout_depsilons_p, breakFlag] = ...
    SPMe_Step_fcn(I,x_cse_n,x_cse_p,x_cen,x_cep,x_dcse_depsilons_p,params)
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
%% SPMe Simulation
dt = params.dt;
%------- ELECTROCHEMICAL COMPUTATION -------% 

%%% Intercalation Current Density
j_Li_n = I/(Area_n*L_n); % [A/m^3] deintercalation current density in negative electrode
j_Li_p = -I/(Area_p*L_p); % [A/m^3] intercalation current density in positive electrode

%%% Electrolyte Concentration (1st-Order Padé Approximation)
% State-Space Representation
a_p0 =  -(861470818036363350116650002266536*epsilon_e_n*epsilon_e_sep^(3/2) + 9408126880630107114529943393067*epsilon_e_n^(3/2)*epsilon_e_sep + 133595401704947513588597985661860*epsilon_e_n^(5/2) + 59271199347969668977610120825136*epsilon_e_sep^(5/2))/(885443715538058477568*D_e*epsilon_e_n^(3/2)*epsilon_e_sep^(3/2)*(23729891576419966*epsilon_e_n + 1770887431076117*epsilon_e_sep));
b_p0 = (10005376170368243727572786501173029165093451132494604498550212087083751263089889024*epsilon_e_n*epsilon_e_sep^4 + 339063137833242250763345947871294659145374370532520699818750855944421000525516000*epsilon_e_n^4*epsilon_e_sep + 1991198357367298664771433826220905659533808478515343956297310707128673204699456000*epsilon_e_n^5 + 326616175708804789147026173844824397788330462258596596588959980085671013139148800*epsilon_e_sep^5 + 70884816238697108610900131788805299059083655277109067516205435967067381596603961547*epsilon_e_n^2*epsilon_e_sep^3 + 17503336390750479458472689125376349464874370860441563495980986971424426872661450*epsilon_e_n^3*epsilon_e_sep^2 + 107421719924885443449939735519594730178119483538388288908706524098701318848852240*epsilon_e_n^(3/2)*epsilon_e_sep^(7/2) + 3124603654649766990761131302908245417082392305115329003801114553555894997984704320*epsilon_e_n^(5/2)*epsilon_e_sep^(5/2) + 21733976812409541213609061507809290224732901039587932414940505083601390030466201600*epsilon_e_n^(7/2)*epsilon_e_sep^(3/2))/(1306684288976403699699358492537989931991040*D_e*epsilon_e_n^(1/2)*epsilon_e_sep^(3/2)*(2932066978031150561007770110743347109988322986088*epsilon_e_n*epsilon_e_sep^(5/2) + 459836248543411033322937842099687545266888573342*epsilon_e_n^(5/2)*epsilon_e_sep + 3170204397566675566700636857099689540137828696760*epsilon_e_n^(7/2) + 104962621950126428042449790681193440832462876912*epsilon_e_sep^(7/2) + 20442609108252715977175690876820142567952604057776*epsilon_e_n^2*epsilon_e_sep^(3/2) + 16660733642877212443077676091160518542127080839*epsilon_e_n^(3/2)*epsilon_e_sep^2));

a_n0 = (827852444649578427085517608056932*epsilon_e_n*epsilon_e_sep^(3/2) + 9408126880630107114529943393067*epsilon_e_n^(3/2)*epsilon_e_sep + 118542398695939337955220241650272*epsilon_e_n^(5/2) + 66797700852473756794298992830930*epsilon_e_sep^(5/2))/(885443715538058477568*D_e*epsilon_e_n^(3/2)*epsilon_e_sep^(3/2)*(23729891576419966*epsilon_e_n + 1770887431076117*epsilon_e_sep));
b_n0 = (11157123710217850623701305410960527293796215865349497556006254291741831436635367456*epsilon_e_n*epsilon_e_sep^4 + 316658867253081630730782458925393769903680827188656997201278631774276871798405760*epsilon_e_n^4*epsilon_e_sep + 1766837979072391697128390764538252514048674654467707605885625949874419409895731200*epsilon_e_n^5 + 467510196339544349427911583802327806974184061324572791295760832064308352849225000*epsilon_e_sep^5 + 70099873283808439733950362070934608353720127973433475288655737045829323602876929003*epsilon_e_n^2*epsilon_e_sep^3 + 17503336390750479458472689125376349464874370860441563495980986971424426872661450*epsilon_e_n^3*epsilon_e_sep^2 + 127958967956699345146456267053337211983005231603596682974722729588000103515369960*epsilon_e_n^(3/2)*epsilon_e_sep^(7/2) + 3161194509996544238763524950835006709420740658279824224696997891997929502236027280*epsilon_e_n^(5/2)*epsilon_e_sep^(5/2) + 20230762277833664644425052479070393761714535556805773599395285728677128700338086400*epsilon_e_n^(7/2)*epsilon_e_sep^(3/2))/(1306684288976403699699358492537989931991040*D_e*epsilon_e_n^(1/2)*epsilon_e_sep^(3/2)*(3051135687798913063854203483921518325063075841424*epsilon_e_n*epsilon_e_sep^(5/2) + 433179074714807485887512055849745023975530529546*epsilon_e_n^(5/2)*epsilon_e_sep + 2812998268263388058161336737565175894913570130752*epsilon_e_n^(7/2) + 118291208864428201760162683806164701478141898810*epsilon_e_sep^(7/2) + 19644848752808707268617659372739256882990069504312*epsilon_e_n^2*epsilon_e_sep^(3/2) + 16660733642877212443077676091160518542127080839*epsilon_e_n^(3/2)*epsilon_e_sep^2));

A_cep = -1/b_p0;
B_cep = 1;
C_cep = a_p0/b_p0;

A_cen = -1/b_n0;
B_cen = 1;
C_cen = a_n0/b_n0;

% [G_cep,H_cep] = c2d(A_cep,B_cep,dt);
% [G_cen,H_cen] = c2d(A_cen,B_cen,dt);
G_cep = 0.9498;
G_cen = 0.9514;
H_cep = 0.9747;
H_cen = 0.9755;
% Update (note that y must be updated before x because we only know the 'current' current value)
y_cep = C_cep*x_cep;
x_cep_updated = G_cep*x_cep + H_cep*gamma*I;
ce_p = y_cep + ce_0;

y_cen = C_cen*x_cen;
x_cen_updated = G_cen*x_cen + H_cen*gamma*I;
ce_n = y_cen + ce_0;

%%% Electrode Surface Concentration (3rd-Order Padé Approximation)
% Negative Electrode: Derived Parameters
m_n = 1/(D_s_n*a_s_n*F); % [s-mol/m-C] convenience variable
a1_cse_n = 189*D_s_n/R_s_n^2;
a2_cse_n = 3465*D_s_n^2/R_s_n^4;
b0_cse_n = -21*m_n*D_s_n/R_s_n;
b1_cse_n = -1260*m_n*D_s_n^2/R_s_n^3;
b2_cse_n = -10395*m_n*D_s_n^3/R_s_n^5;
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
A_cse_p = [0 1 0; 0 0 1; 0 -a2_cse_p -a1_cse_p]; % A matrix
B_cse_p = [0 0 1]'; % B vector
C_cse_p = [b2_cse_p b1_cse_p b0_cse_p]; % C vector

% [G_cse_p,H_cse_p] = c2d(A_cse_p,B_cse_p,dt);
% [G_cse_n,H_cse_n] = c2d(A_cse_n,B_cse_n,dt);
G_cse_p = [1,0.999987251922180,0.495333120125527;0,0.999961845294749,0.986028939328108;0,-7.59522067480488e-05,0.972175897610784];
H_cse_p = [0.165498465272271;0.495333120125527;0.986028939328108];
G_cse_n = [1,0.897892611033907,0.184383559640740;0,0.758107772241298,0.219810574537802;0,-0.288368820210808,-0.0502593203367288];
H_cse_n = [0.0778318675950616;0.184383559640740;0.219810574537802];
% Update (note that y must be updated before x because we only know the 'current' current value)
cse_n = C_cse_n*x_cse_n; % [mol/m^3] negative electrode surface concentration
x_cse_n_updated = G_cse_n*x_cse_n + H_cse_n*j_Li_n; % update state variables
cse_p = C_cse_p*x_cse_p; % [mol/m^3] positive electrode surface concentration
x_cse_p_updated = G_cse_p*x_cse_p + H_cse_p*j_Li_p; % update state variables

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
V_out = U_p - U_n + dphi_econ + eta_p - eta_n - I*R_lump; % [V] terminal voltage

%%% SOC
% SOC_n = (theta_n - theta_n_0)/(theta_n_1 - theta_n_0); % [-] SOC from negative electrode surface stoichiometry
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
% Positive Electrode: State-Space Representation
% A_dcse_depsilons_p = A_cse_p; % A matrix
% B_dcse_depsilons_p = B_cse_p; % B vector
C_dcse_depsilons_p = -1/epsilon_s_p*C_cse_p; % C vector

% [G_dcse_depsilons_p,H_dcse_depsilons_p] = c2d(A_dcse_depsilons_p,B_dcse_depsilons_p,dt);
G_dcse_depsilons_p = [1,0.999987251922180,0.495333120125527;0,0.999961845294749,0.986028939328108;0,-7.59522067480488e-05,0.972175897610784];
H_dcse_depsilons_p = [0.165498465272271;0.495333120125527;0.986028939328108];
% Update
dcse_depsilons_p = C_dcse_depsilons_p*x_dcse_depsilons_p; % [mol/m^3] derivative of positive electrode surface concentration wrt active material volume fraction 
x_dcse_depsilons_p_updated = G_dcse_depsilons_p*x_dcse_depsilons_p + H_dcse_depsilons_p*j_Li_p; % update state variables


%%% Terminal Voltage Sensitivity WRT Active Material Volume Fraction
dVout_depsilons_p = -dRlump_depsilons_p*I+deta_das_das_depsilons_p+(deta_dcse_p+dU_dcse_p)*dcse_depsilons_p; % [V] derivative of terminal voltage wrt positive electrode active material volume fraction


%------- STOP CONDITION -------% 

cse_p_next = C_cse_p*x_cse_p_updated; % check next cse_p for problems
cse_n_next = C_cse_n*x_cse_n_updated; % check next cse_n for problems
breakFlag = 0; % initialize break flag
if V_out < V_out_min || V_out > V_out_max || SOC < SOC_min || SOC > SOC_max || cse_p_next > cse_p_max || cse_p_next < 0 || cse_n_next > cse_n_max || cse_n_next < 0 % if voltage or SOC or ionic concentration limit is reached
    breakFlag = 1; % raise flag to indicate that the loop that is iterating this function should break
end

end