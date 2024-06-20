%%%--------Implementation of Q-learning--------%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc
rng('default') % reset random number generator

%% Read the .ini file
ini = IniConfig();
ini.ReadFile('parameters.ini');
sections = ini.GetSections();

%% User-defined parameters
% get values from .ini files 
[keys, count_keys] = ini.GetKeys(sections{1});
values = ini.GetValues(sections{1}, keys);
% battery electrochemical parameters
params.epsilon_s_n = values{1}; % [-] active material volume fraction of solid negative electrode 
params.epsilon_s_p = values{2}; % [-] active material volume fraction of solid positive electrode 
params.theta_n_0 = values{3}; % [-] negative electrode surface stoichiometry when SOC = 0 
params.theta_n_1 = values{4}; % [-] negative electrode surface stoichiometry when SOC = 1 
params.theta_p_0 = values{5}; % [-] positive electrode surface stoichiometry when SOC = 0 
params.theta_p_1 = values{6}; % [-] positive electrode surface stoichiometry when SOC = 1   
params.D_s_n = str2num(values{7}); % [m^2/s] solid-phase diffusion coefficient for negative electrode
params.D_s_p = str2num(values{8}); % [m^2/s] solid-phase diffusion coefficient for positive electrode
params.R_s_n = str2num(values{9}); % [m] radius of solid negative electrode particle
params.R_s_p = str2num(values{10}); % [m] radius of solid positive electrode particle
params.F = values{11}; % [C/mol] Faraday constant
params.cse_n_max = values{12}; % [mol/m^3] max concentration in negative electrode
params.cse_p_max = values{13}; % [mol/m^3] max concentration in positive electrode
params.epsilon_e_sep = values{14}; % [-] active material volume fraction (porosity) of electrolyte in separator
params.epsilon_e_n = values{15}; % [-] active material volume fraction (porosity) of electrolyte in negative electrode
params.epsilon_e = values{16}; % [-] uniform active material volume fraction of electrolyte in electrodes (for Pade approximation)
params.k_n = str2num(values{17}); % [(A/m^2)*(m^3/mol)^(1-1/(3*alpha))] reaction rate in negative electrode
params.k_p = str2num(values{18}); % [(A/m^2)*(m^3/mol)^(1-1/(3*alpha))] reaction rate in positive electrode
params.kappa = values{19};  % [1/Ohms-m] electrolyte conductivity
params.t_plus = values{20}; % [-] transference number
params.D_e = str2num(values{21}); % [m^2/s] electrolyte diffusion coefficient (assumed constant for all phases and x-coordinates; from Capiglia et al., 1999)
params.L_n = str2num(values{22}); % [m] thickness of negative electrode
params.L_sep = str2num(values{23}); % [m] thickness of separator
params.L_p = str2num(values{24}); % [m] thickness of positive electrode
params.Area_n = values{25}; % [m^2] negative electrode current collector area
params.Area_p = values{26}; % [m^2] positive electrode current collector area
params.R = values{27}; % [J/mol-K] universal gas constant
params.T = values{28}; % [K] cell temperature
params.alpha = values{29}; % [-] charge transfer coefficient (for positive and negative electrodes)
params.R_f_n = values{30}; % [Ohms*m^2] resistivity of SEI layer for negative electrode
params.R_f_p = values{31}; % [Ohms*m^2] resistivity of SEI layer for positive electrode
params.R_c = values{32}; %  [Ohms*m^2] contact resistance (resistance of current collectors and wiring)

% get values from .ini files
[keys, count_keys] = ini.GetKeys(sections{2});
values = ini.GetValues(sections{2}, keys);
% Limits
params.V_out_max = values{1}; % [V] terminal voltage upper limit
params.V_out_min = values{2}; % [V] terminal voltage lower limit
params.SOC_max = values{3}; % [-] SOC upper limit
params.SOC_min = values{4}; % [-] SOC lower limit
params.I_max = values{5};   % [A] Current input upper limit
params.I_min = values{6};  % [A] Current input lower limit
params.dcse_depsilonsp_max = values{7};    % [mol/m^3] upper limit of the normalized state sensitivity dcse/depsilon_sp
params.dcse_depsilonsp_min = values{8};   % [mol/m^3] lower limit of the normalized state sensitivity dcse/depsilon_sp

% get values from .ini files
[keys, count_keys] = ini.GetKeys(sections{3});
values = ini.GetValues(sections{3}, keys);
% initial battery states
params.SOC_0 = values{1}; % [-] battery initial SOC

% get values from .ini files
[keys, count_keys] = ini.GetKeys(sections{4});
values = ini.GetValues(sections{4}, keys);
% Q-learning initialization
params.t_final = values{1};           % optimization time horizon
params.dt = values{2};                   % length of time step 
params.n_S = values{3};                % number of agent states discretization
params.n_action = values{4};            % number of action discretization
params.n_episode = values{5};      % number of episodes
params.epsilon_ini = values{6};          % greedy parameter epsilon at the beginning of the episode
params.epsilon_end = values{7};          % greedy parameter epsilon at the end of the episode
params.alpha_ini = values{8};          % learning rate alpha at the beginning of the episode
params.alpha_end = values{9};         % learning rate alpha at the beginning of the episode
params.gamma = values{10};            % decaying factor gamma

%% Get file directory values
[keys, count_keys] = ini.GetKeys(sections{5});
values = ini.GetValues(sections{5}, keys);
root_directory = convertCharsToStrings(values{1});

%% Computed parameters
params.a_s_n = 3*params.epsilon_s_n/params.R_s_n; % [1/m] negative electrode specific interfacial surface area
params.a_s_p = 3*params.epsilon_s_p/params.R_s_p; % [1/m] positive electrode specific interfacial surface area

params.m_p = 1/(params.D_s_p*params.a_s_p*params.F);     % [s-mol/m-C] convenience coefficient for positive particle
params.b0_cse_p = -21*params.m_p*params.D_s_p/params.R_s_p;        
params.b1_cse_p = -1260*params.m_p*params.D_s_p^2/params.R_s_p^3;
params.b2_cse_p = -10395*params.m_p*params.D_s_p^3/params.R_s_p^5;
params.C_cse_p = [params.b2_cse_p params.b1_cse_p params.b0_cse_p];      % output matrix C for cse_p
params.C_dcse_depsilons_p = -1/params.epsilon_s_p*params.C_cse_p; % output matrix C for dcse_p/depsilon_sp

params.m_n = 1/(params.D_s_n*params.a_s_n*params.F); % [s-mol/m-C] convenience coefficient for positive particle
params.b2_cse_n = -10395*params.m_n*params.D_s_n^3/params.R_s_n^5;   

%% battery model initialization
theta_n = (params.theta_n_1-params.theta_n_0)*params.SOC_0+params.theta_n_0; % [-] initial negative electrode surface stoichiometry
theta_p = params.theta_p_0-(params.theta_p_0-params.theta_p_1)*params.SOC_0; % [-] initial positive electrode surface stoichiometry
cse_0_n = theta_n*params.cse_n_max; % [mol/m^3] initial concentration of ions in negative solid electrode
cse_0_p = theta_p*params.cse_p_max; % [mol/m^3] initial concentration of ions in positive solid electrode

%% Q-learning initialization
time = 0:params.dt:params.t_final;      % time series
Q = zeros(params.n_S,params.n_action);  % initial value of state-action pair
epCount = 0;              % initializing the counter of episodes  
G_max = 0;                % initializing the maximum return

%% discretized actions
actions = params.I_min:(params.I_max-params.I_min)/(params.n_action-1):params.I_max;        % series of discretized action       
S_S = params.dcse_depsilonsp_min:(params.dcse_depsilonsp_max-params.dcse_depsilonsp_min)/(params.n_S-1):params.dcse_depsilonsp_max; % series of discretized agent atate

%% Q-learning: loop for each episode
for i = 1: params.n_episode
    G = 0;             % initialize the 
    SOC = params.SOC_0;       % initialize starting SOC
    time_cnt = 1;      % initialize the counter of steps in an episode    
    states_ep = [];
    actions_ep = [];
    
    x_cse_n = zeros(3,1); % negative electrode: initialize state variables
    x_cse_n(1) = cse_0_n/params.b2_cse_n; % negative electrode: initial condition for state variables (uniform concentration profile)
    x_cse_p = zeros(3,1); % positive electrode: initialize state variables
    x_cse_p(1) = cse_0_p/params.b2_cse_p; % positive electrode: initial condition for state variables (uniform concentration profile)
    x_cen = 0; % state variable initial condition
    x_cep = 0;
    x_dcse_depsilons_p = zeros(3,1); % positive electrode: initialize state variables
    
    dcse_depsilons_p = params.C_dcse_depsilons_p*x_dcse_depsilons_p*params.epsilon_s_p;
    [~,S0_index] = min(abs(dcse_depsilons_p-S_S));
    S = S0_index;
      
    alpha = (params.alpha_ini - params.alpha_end)*(params.n_episode-i)/params.n_episode + params.alpha_end;  % linearly decaying alpha
    epsilon = (params.epsilon_ini - params.epsilon_end)*(params.n_episode-i)/params.n_episode + params.epsilon_end;  % linearly decaying
    
    % loop for each step of episode 
    while time_cnt < length(time)
        
        % choosing A from S using policy derived from Q (epsilon-greedy)
        if rand < epsilon
            A = actions(randi(length(actions)));
        else
            A = actions(Q(S,:) == max(Q(S,:)));
            A = A(randi(length(A)));
        end
        
        states_ep = [states_ep SOC];
        actions_ep = [actions_ep A];
        x_ce = 1;
        % Take action A, observe R S'        
        % update battery states
        [x_cse_n_updated, x_cse_p_updated, x_cen_updated, x_cep_updated, x_dcse_depsilons_p_updated, V_out,SOC,dVout_depsilons_p, breakFlag]...
            = SPMe_Step_fcn(A,x_cse_n,x_cse_p,x_cen,x_cep,x_dcse_depsilons_p,params);
        x_cse_n = x_cse_n_updated;
        x_cse_p = x_cse_p_updated;
        x_cen = x_cen_updated;
        x_cep = x_cep_updated;
        x_dcse_depsilons_p = x_dcse_depsilons_p_updated;
        dcse_depsilons_p = params.C_dcse_depsilons_p*x_dcse_depsilons_p*params.epsilon_s_p;
            
        if breakFlag == 1            % if the voltage get out of the range
            R = -1;
            Q(S,A==actions) = Q(S,A==actions) + alpha*(R+params.gamma*0 - Q(S,A==actions));
            break

        else
            R = (dVout_depsilons_p*params.epsilon_s_p)^2; % observe intermediate reward
            time_cnt_next = time_cnt+1; % increment state
            [~,S_index] = min(abs(dcse_depsilons_p-S_S)); % discretize V_out (find index of nearest value in S_Vout)
            S_next = S_index; % next terminal voltage state
            % Q(S,A) <-- Q(S,A)+alpha[R+gamma*max{Q(S',A')}-Q(S,A)]
            Q(S,A==actions) = Q(S,A==actions) + alpha*(R+params.gamma*max(Q(S_next,:)) - Q(S,A==actions));
        end
        
        G = G + R;
        time_cnt = time_cnt_next;
        S = S_next; 
    end
    
    epCount = epCount+1; % increment episode counter
    G_log(epCount) = G;                    
    
    if G > G_max
        G_max = G;
        states_max = states_ep;
        actions_max = actions_ep;
    end     
end

%% Generate learned policy
pi = [];     % initialize policy
i = 1;       % initialize counter
G = 0;       % initialize return for episode

% Initial state:
s_time = 1;                    % initial state is first time step
x_cse_n = zeros(3,1);          % negative electrode: initialize state variables
x_cse_n(1) = cse_0_n/params.b2_cse_n; % negative electrode: initial condition for state variables (uniform concentration profile)
x_cse_p = zeros(3,1);          % positive electrode: initialize state variables
x_cse_p(1) = cse_0_p/params.b2_cse_p; % positive electrode: initial condition for state variables (uniform concentration profile)
x_cen = 0;                     % state variable initial condition
x_cep = 0;
x_dcse_depsilons_p = zeros(3,1); % positive electrode: initialize state variables
dcse_depsilons_p = params.C_dcse_depsilons_p*x_dcse_depsilons_p*params.epsilon_s_p;

[~,S0_index] = min(abs(dcse_depsilons_p-S_S));
S = S0_index;

pi(i) = 0;
SOC_plot(i) = params.SOC_0;

while s_time < length(time) % for each nonterminal state of the episode
    i = i+1; % increment counter
    
    a = actions(Q(S,:) == max(Q(S,:))); % determine action(s) that maximize Q
    if length(a) > 1 % if there are multiple actions that maximize Q
        a = a(randi(length(a))); % break ties randomly
    end

    % Take action a, observe R and s_next:
    [x_cse_n_updated, x_cse_p_updated, x_cen_updated, x_cep_updated, x_dcse_depsilons_p_updated, V_out,SOC, dVout_depsilons_p, breakFlag] = ...
        SPMe_Step_fcn(a,x_cse_n,x_cse_p,x_cen, x_cep, x_dcse_depsilons_p,params);
    x_cse_n = x_cse_n_updated; % update states
    x_cse_p = x_cse_p_updated;
    x_cen = x_cen_updated;
    x_cep = x_cep_updated;
    x_dcse_depsilons_p = x_dcse_depsilons_p_updated;
    G = G + (dVout_depsilons_p*params.epsilon_s_p)^2; % update return (sum of normalized sensitivities squared)
    
    dcse_depsilons_p = params.C_dcse_depsilons_p*x_dcse_depsilons_p*params.epsilon_s_p;
    % Update pi:
    pi(i) = a; % update pi
    SOC_plot(i) = SOC;
    % Update s:
    if breakFlag % if action is invalid
        disp("Policy is invalid")
        return
    else % intermediate state with valid action
        s_time = s_time+1; % increment state   
        [~,S_index] = min(abs(dcse_depsilons_p-S_S)); % discretize V_out (find index of nearest value in S_Vout)
        S = S_index; % next terminal voltage state
    end
end

%% Save the active learning results
I = pi;
writematrix(Q,root_directory+"\Q_table.csv");
writematrix([time' I'],root_directory+"\learned_policy.csv");
