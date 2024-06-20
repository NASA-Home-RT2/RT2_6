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
params.w2_max = values{7}; % [V] Upper bound for the augmented state w2
params.w2_min = values{8}; % [V] Uower bound for the augmented state w2
params.Time_max = values{9};
params.Time_min = values{10};

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
params.n_time = values{11};
params.n_w2 = values{12};

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

params.a_p0 =  -(861470818036363350116650002266536*params.epsilon_e_n*params.epsilon_e_sep^(3/2) + 9408126880630107114529943393067*params.epsilon_e_n^(3/2)*params.epsilon_e_sep + 133595401704947513588597985661860*params.epsilon_e_n^(5/2) + 59271199347969668977610120825136*params.epsilon_e_sep^(5/2))/(885443715538058477568*params.D_e*params.epsilon_e_n^(3/2)*params.epsilon_e_sep^(3/2)*(23729891576419966*params.epsilon_e_n + 1770887431076117*params.epsilon_e_sep));
params.b_p0 = (10005376170368243727572786501173029165093451132494604498550212087083751263089889024*params.epsilon_e_n*params.epsilon_e_sep^4 + 339063137833242250763345947871294659145374370532520699818750855944421000525516000*params.epsilon_e_n^4*params.epsilon_e_sep + 1991198357367298664771433826220905659533808478515343956297310707128673204699456000*params.epsilon_e_n^5 + 326616175708804789147026173844824397788330462258596596588959980085671013139148800*params.epsilon_e_sep^5 + 70884816238697108610900131788805299059083655277109067516205435967067381596603961547*params.epsilon_e_n^2*params.epsilon_e_sep^3 + 17503336390750479458472689125376349464874370860441563495980986971424426872661450*params.epsilon_e_n^3*params.epsilon_e_sep^2 + 107421719924885443449939735519594730178119483538388288908706524098701318848852240*params.epsilon_e_n^(3/2)*params.epsilon_e_sep^(7/2) + 3124603654649766990761131302908245417082392305115329003801114553555894997984704320*params.epsilon_e_n^(5/2)*params.epsilon_e_sep^(5/2) + 21733976812409541213609061507809290224732901039587932414940505083601390030466201600*params.epsilon_e_n^(7/2)*params.epsilon_e_sep^(3/2))/(1306684288976403699699358492537989931991040*params.D_e*params.epsilon_e_n^(1/2)*params.epsilon_e_sep^(3/2)*(2932066978031150561007770110743347109988322986088*params.epsilon_e_n*params.epsilon_e_sep^(5/2) + 459836248543411033322937842099687545266888573342*params.epsilon_e_n^(5/2)*params.epsilon_e_sep + 3170204397566675566700636857099689540137828696760*params.epsilon_e_n^(7/2) + 104962621950126428042449790681193440832462876912*params.epsilon_e_sep^(7/2) + 20442609108252715977175690876820142567952604057776*params.epsilon_e_n^2*params.epsilon_e_sep^(3/2) + 16660733642877212443077676091160518542127080839*params.epsilon_e_n^(3/2)*params.epsilon_e_sep^2));

params.a_n0 = (827852444649578427085517608056932*params.epsilon_e_n*params.epsilon_e_sep^(3/2) + 9408126880630107114529943393067*params.epsilon_e_n^(3/2)*params.epsilon_e_sep + 118542398695939337955220241650272*params.epsilon_e_n^(5/2) + 66797700852473756794298992830930*params.epsilon_e_sep^(5/2))/(885443715538058477568*params.D_e*params.epsilon_e_n^(3/2)*params.epsilon_e_sep^(3/2)*(23729891576419966*params.epsilon_e_n + 1770887431076117*params.epsilon_e_sep));
params.b_n0 = (11157123710217850623701305410960527293796215865349497556006254291741831436635367456*params.epsilon_e_n*params.epsilon_e_sep^4 + 316658867253081630730782458925393769903680827188656997201278631774276871798405760*params.epsilon_e_n^4*params.epsilon_e_sep + 1766837979072391697128390764538252514048674654467707605885625949874419409895731200*params.epsilon_e_n^5 + 467510196339544349427911583802327806974184061324572791295760832064308352849225000*params.epsilon_e_sep^5 + 70099873283808439733950362070934608353720127973433475288655737045829323602876929003*params.epsilon_e_n^2*params.epsilon_e_sep^3 + 17503336390750479458472689125376349464874370860441563495980986971424426872661450*params.epsilon_e_n^3*params.epsilon_e_sep^2 + 127958967956699345146456267053337211983005231603596682974722729588000103515369960*params.epsilon_e_n^(3/2)*params.epsilon_e_sep^(7/2) + 3161194509996544238763524950835006709420740658279824224696997891997929502236027280*params.epsilon_e_n^(5/2)*params.epsilon_e_sep^(5/2) + 20230762277833664644425052479070393761714535556805773599395285728677128700338086400*params.epsilon_e_n^(7/2)*params.epsilon_e_sep^(3/2))/(1306684288976403699699358492537989931991040*params.D_e*params.epsilon_e_n^(1/2)*params.epsilon_e_sep^(3/2)*(3051135687798913063854203483921518325063075841424*params.epsilon_e_n*params.epsilon_e_sep^(5/2) + 433179074714807485887512055849745023975530529546*params.epsilon_e_n^(5/2)*params.epsilon_e_sep + 2812998268263388058161336737565175894913570130752*params.epsilon_e_n^(7/2) + 118291208864428201760162683806164701478141898810*params.epsilon_e_sep^(7/2) + 19644848752808707268617659372739256882990069504312*params.epsilon_e_n^2*params.epsilon_e_sep^(3/2) + 16660733642877212443077676091160518542127080839*params.epsilon_e_n^(3/2)*params.epsilon_e_sep^2));

params.C_cep = params.a_p0/params.b_p0*(params.Area_p*params.L_p*params.F*params.epsilon_s_p);
params.C_cen = params.a_n0/params.b_n0*(params.Area_n*params.L_n*params.F*params.epsilon_s_n);

params.A_cep = -1/params.b_p0*(params.R_s_p^2/params.D_s_p);
params.B_cep = 1;

params.A_cen = -1/params.b_n0*(params.R_s_n^2/params.D_s_n);
params.B_cen = 1;

params.dt_non_dim = params.D_s_p / params.R_s_p^2;
[params.G_cep,params.H_cep] = c2d(params.A_cep,params.B_cep,params.dt_non_dim);
[params.G_cen,params.H_cen] = c2d(params.A_cen,params.B_cen,params.dt_non_dim);

params.A_cse_p = [0 1 0 ; 0 0 1; 0 -3465 -189];
params.B_cse_p = [0; 0; 1];
params.A_dcse_dDs_p = [0 1 0 0; 0 0 1 0; 0 0 0 1; -12006225 -1309770 -42651 -378];
params.B_dcse_dDs_p = [0 0 0 1]';
[params.G_cse_p,params.H_cse_p] = c2d(params.A_cse_p,params.B_cse_p,params.dt_non_dim);
[params.G_dcse_dDs_p,params.H_dcse_dDs_p] = c2d(params.A_dcse_dDs_p,params.B_dcse_dDs_p,params.dt_non_dim);

%% Non-dimensional scaling factors
scale_I = params.R_s_p.^2/params.D_s_p/params.F/params.epsilon_s_p/params.Area_p/params.L_p; % [-] Sacling factor to convert I to \tilde I

%% battery model initialization
theta_n = (params.theta_n_1-params.theta_n_0)*params.SOC_0+params.theta_n_0; % [-] initial negative electrode surface stoichiometry
theta_p = params.theta_p_0-(params.theta_p_0-params.theta_p_1)*params.SOC_0; % [-] initial positive electrode surface stoichiometry
cse_0_n = theta_n*params.cse_n_max; % [mol/m^3] initial concentration of ions in negative solid electrode
cse_0_p = theta_p*params.cse_p_max; % [mol/m^3] initial concentration of ions in positive solid electrode

%% Q-learning initialization
Q = zeros(params.n_S,params.n_w2,params.n_time,params.n_action);  % initial value of state-action pair
epCount = 0;              % initializing the counter of episodes  
G_max = 0;                % initializing the maximum return

%% discretized actions and states
actions = (params.I_min:(params.I_max-params.I_min)/(params.n_action-1):params.I_max)*scale_I;        % series of discretized action       
States_SOC = params.SOC_min:(params.SOC_max-params.SOC_min)/(params.n_S-1):params.SOC_max;     % series of discretized agent state SOC
States_w2 = params.w2_min:(params.w2_max-params.w2_min)/(params.n_w2-1):params.w2_max;         % series of discretized agent state w2

Dt = params.t_final/params.n_time;    % 
time = 0:Dt:params.t_final;           % time series

%% Q-learning: loop for each episode
for i = 1: params.n_episode

    G = 0;                    % initialize the return
    SOC = params.SOC_0;       % initialize starting SOC
    w2 = 0;                   % initialize starting w2
    R1 = 0; 
    R2 = 0;

    time_cnt = 1;      % initialize the counter of steps in an episode    
    states_ep = [];
    actions_ep = [];
   
    x_cse_p = zeros(3,1); % positive electrode: initialize state variables
    x_cse_p(1) = cse_0_p/3465; % positive electrode: initial condition for state variables (uniform concentration profile)
    x_cen = 0; % state variable initial condition
    x_cep = 0;
    x_dcse_dDs_p = zeros(4,1); % positive electrode: initialize state variables
    
    
    [~,S_SOC] = min(abs(SOC-States_SOC));
    [~,S_w2] = min(abs(w2-States_w2));

    alpha = params.alpha_ini;  % fixed alpha
    epsilon = (params.epsilon_ini - params.epsilon_end)*(params.n_episode-i)/params.n_episode + params.epsilon_end;  % linearly decaying
    
    % loop for each step of episode 
    while time_cnt < length(time)
        
        % choosing A from S using policy derived from Q (epsilon-greedy)
        if rand < epsilon
            A = actions(randi(length(actions)));
        else
            A = actions(Q(S_SOC,S_w2,time_cnt,:) == max(Q(S_SOC,S_w2,time_cnt,:)));
            A = A(randi(length(A)));
        end
        
        states_ep = [states_ep SOC];
        actions_ep = [actions_ep A];
      
        % Take action in Dt steps        
        % update battery states
        [x_cse_p_updated, x_cen_updated, x_cep_updated, x_dcse_dDs_p_updated, V_out,SOC,dVout_depsilons_p,dw2,dR1,dR2, breakFlag]...
            = SPMe_Step_fcn(A,x_cse_p,x_cen,x_cep,x_dcse_dDs_p,Dt,params);
        x_cse_p = x_cse_p_updated;
        x_cen = x_cen_updated;
        x_cep = x_cep_updated;
        x_dcse_dDs_p = x_dcse_dDs_p_updated;
                    
        if breakFlag == 1            % if the voltage get out of the range
            w2 = w2 + dw2;
            R_diff = (R1 + dR1)*(R2 + dR2) - R1*R2;
            R = R_diff - w2^2 - 1e9;    
            Q(S_SOC,S_w2,time_cnt,A==actions) = R;
            G = G + R;
            break
        else
            w2 = w2 + dw2;
            R_diff = (R1 + dR1)*(R2 + dR2) - R1*R2;
            R1 = R1 + dR1;
            R2 = R2 + dR2;
            % update state indexes
            [~,S_SOC_next] = min(abs(SOC-States_SOC));        
            [~,S_w2_next] = min(abs(w2-States_w2));
            time_cnt_next = time_cnt+1;        

            R = R_diff - floor(time_cnt/params.n_time)*w2^2;
            % Q(S,A) <-- Q(S,A)+alpha[R+gamma*max{Q(S',A')}-Q(S,A)]
            if time_cnt < params.n_time
                Q(S_SOC,S_w2,time_cnt,A==actions) = Q(S_SOC,S_w2,time_cnt,A==actions) + ...
                    alpha*( R + params.gamma*max(Q(S_SOC_next,S_w2_next,time_cnt_next,:)) - Q(S_SOC,S_w2,time_cnt,A==actions));
            else
                Q(S_SOC,S_w2,time_cnt,A==actions) = R;
            end
        end
        
        G = G + R;
        time_cnt = time_cnt_next;
        S_SOC = S_SOC_next;
        S_w2 = S_w2_next;
    end
    
    epCount = epCount+1; % increment episode counter             
    
end

%% Generate learned policy
t_final_index = params.t_final;
Dt = t_final_index/(params.n_time);
time = 0:Dt:t_final_index; % time series

I = [];
SOC_plot = [];
s_time = 1; % initial state is first time step

% Initial state:
SOC_0 = params.SOC_0;
R1 = 0;
R2 = 0;
w2 = 0;

[~,S_SOC] = min(abs(SOC_0-States_SOC));
[~,S_w2] = min(abs(w2-States_w2));

x_cse_p = zeros(3,1);          % positive electrode: initialize state variables
x_cse_p(1) = cse_0_p/3465;     % positive electrode: initial condition for state variables (uniform concentration profile)
x_cen = 0;                     % state variable initial condition
x_cep = 0;
x_dcse_dDs_p = zeros(4,1); % positive electrode: initialize state variables

SOC_plot(s_time) = SOC_0;

while s_time < length(time) % for each nonterminal state of the episode
    
    a = actions(Q(S_SOC,S_w2,s_time,:) == max(Q(S_SOC,S_w2,s_time,:))); % determine action(s) that maximize Q
    if length(a) > 1 % if there are multiple actions that maximize Q
        a = a(randi(length(a))); % break ties randomly
    end

    % Take action a, observe R and s_next:
    [x_cse_p_updated, x_cen_updated, x_cep_updated, x_dcse_dDs_p_updated, V_out,SOC,dVout_depsilons_p,dw2,dR1,dR2, breakFlag]...
            = SPMe_Step_fcn(A,x_cse_p,x_cen,x_cep,x_dcse_dDs_p,Dt,params);
    
    % update diffusion states
    x_cse_p = x_cse_p_updated; % update states
    x_cen = x_cen_updated;
    x_cep = x_cep_updated;
    x_dcse_dDs_p = x_dcse_dDs_p_updated;

    R1 = R1 + dR1;
    R2 = R2 + dR2;
    w2 = w2 + dw2;

    % Update pi:
    I(s_time) = a; % update pi
    SOC_plot(s_time) = SOC; 
    
    % Update s:
    if breakFlag % if action is invalid
        disp("Policy is invalid")
        return
    else % intermediate state with valid action
        s_time = s_time+1; % increment state   
        [~,S_SOC] = min(abs(SOC-States_SOC)); % discretize V_out (find index of nearest value in S_Vout)
        [~,S_w2] = min(abs(w2-States_w2));    
    end
end

% Plots and generate the time-series of 600s
I_plot = [];
time_plot = 1:params.t_final;
for j = 1:length(I)    
    I_plot = [I_plot I(j)*ones(1,Dt)];
end
%% Save the active learning results
I = I_plot;
writematrix(Q,root_directory+"\Q_table.csv");
writematrix([time_plot' I_plot'],root_directory+"\learned_policy.csv");
