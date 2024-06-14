%% Data Selection Parameter Estimation for SPMe Battery Model (Target Parameter: D_s_p)
% Jack Fogelquist
% 2024-06-04
clc; close all; clear all;

% This script is an example implementation of the data selection parameter
% estimation framework developed in J. Fogelquist and X. Lin,
% "Uncertainty-Aware Data Selection Framework for Parameter Estimation with
% Application to Li-ion Battery," 2022 American Control Conference (ACC),
% Atlanta, GA, USA, 2022, pp. 384-391, doi: 10.23919/ACC53348.2022.9867243

% In this implementation, the framework is applied to estimate a
% health-related electrochemical battery parameter called the cathode
% diffusion coefficient (D_s_p). The estimation voltage data is simulated
% from a random operational input current profile called the Federal Urban
% Driving Schedule (FUDS).

% Algorithm Summary: 
% This algorithm loads a random input current profile; simulates the
% battery output voltage under the true parameter set (for reference) and
% an uncertain parameter set (for estimation); rates the quality of the
% data segments a priori (i.e., under the initial guess of the target
% parameter); extracts high quality data segments and uses them to perform
% the estimation and compute an a posteriori quality rating (i.e., under
% each estimate of the target parameter); and outputs the final estimate
% (based on the a posteriori rating).

% Algorithm Steps
% 1) Load input current profile
% 2) Run SPMe under true parameter values
% 3) Run SPMe under uncertain parameter values
% 4) Build a priori data quality rating table (compute a priori rating for every data segment)
% 5) Build squared errors table (compute squared voltage errors under different values of the target parameter)
% 6) Sweep through the data segments with the lowest a priori ratings and perform the estimation, run the SPMe under the estimate, and compute the a posteriori rating
% 7) Perform the final estimation with the a posteriori data quality rating
% 8) Plot the results

%% Inputs
% Read Inputs from Config File
ini = IniConfig(); % create ini object
ini.ReadFile('config.ini'); % read config file
sections = ini.GetSections(); % load section names
keys = ini.GetKeys(sections{1}); % load key names
values = ini.GetValues(sections{1}, keys); % load key values

% Assign Inputs from Config File
profile_filename = values{1}; % name of input current profile time series data
loadQtable = values{2}; % toggle to load Q-Table [0 = no; 1 = yes]
loadDspErrorTable = values{3}; % toggle to load Dsp Error Table [0 = no; 1 = yes]
numDsp = values{4}; % number of Dsp values to run in SPMe for estimation
uncertainty_Dsp = values{5}; % [%] uncertainty in Dsp
uncertainty_epsilonsp = values{6}; % [%] uncertainty in epsilon_s_p
SOC_0 = values{7}; % [-] initial SOC for battery model

%% Generate Estimation Data
% Load Input Current Profile
profile = load(profile_filename); % load current profile
time = profile.time; % [s] time vector
I = profile.I; % [A] current vector
dt = time(2)-time(1); % [s] time increment

% Model Parameters
R_s_n = 10e-6; % [m] radius of solid negative electrode particle
R_s_p = 10e-6; % [m] radius of solid positive electrode particle
epsilon_s_n = 0.6; % [-] active material volume fraction of solid negative electrode
epsilon_s_p_true = 0.5; % [-] active material volume fraction of solid positive electrode
F = 96485.33289; % [C/mol] Faraday constant
D_s_n = 3.9e-14; % [m^2/s] solid-phase diffusion coefficient for negative electrode
D_s_p_true = 1e-13; % [m^2/s] solid-phase diffusion coefficient for positive electrode

% Computed Parameters
a_s_n = 3*epsilon_s_n/R_s_n; % [1/m] negative electrode specific interfacial surface area
a_s_p = 3*epsilon_s_p_true/R_s_p; % [1/m] positive electrode specific interfacial surface area

% Limits
cse_n_max = 3.6e3 * 372 * 1800 / F; % [mol/m^3] max concentration in negative electrode
cse_p_max = 3.6e3 * 274 * 5010 / F; % [mol/m^3] max concentration in positive electrode
theta_n_0 = 0.0069; % [-] negative electrode surface stoichiometry when SOC = 0
theta_n_1 = 0.6760; % [-] negative electrode surface stoichiometry when SOC = 1
theta_p_0 = 0.8228; % [-] positive electrode surface stoichiometry when SOC = 0
theta_p_1 = 0.4420; % [-] positive electrode surface stoichiometry when SOC = 1 
V_out_max = 4.15; % [V] terminal voltage upper limit
V_out_min = 2.75; % [V] terminal voltage lower limit
SOC_max = 1; % [-] SOC upper limit
SOC_min = 0; % [-] SOC lower limit

% Initial Conditions
ce_0 = 1000; % [mol/m^3] initial concentration of Li-ions in electrolyte
theta_n(1) = (theta_n_1-theta_n_0)*SOC_0+theta_n_0; % [-] initial negative electrode surface stoichiometry
theta_p(1) = theta_p_0-(theta_p_0-theta_p_1)*SOC_0; % [-] initial positive electrode surface stoichiometry
cse_0_n = theta_n(1)*cse_n_max; % [mol/m^3] initial concentration of ions in negative solid electrode
cse_0_p = theta_p(1)*cse_p_max; % [mol/m^3] initial concentration of ions in positive solid electrode

%%% Run SPMe under True Parameters
% Initialize
x_cse_n_step = [cse_0_n/(-10395*D_s_n^2/(a_s_n*F*R_s_n^5)); 0; 0]; % negative electrode: initial condition for state variables (uniform concentration profile)
x_cse_p_step = [cse_0_p/(-10395*D_s_p_true^2/(a_s_p*F*R_s_p^5)); 0; 0]; % positive electrode: initial condition for state variables (uniform concentration profile)
x_ce_step = 0; % state variable initial condition
x_dcse_dDs_n = zeros(4,1); % negative electrode: initialize state variables
x_dcse_dDs_p = zeros(4,1); % positive electrode: initialize state variables
x_dcse_depsilons_n = zeros(3,1); % negative electrode: initialize state variables
x_dcse_depsilons_p = zeros(3,1); % positive electrode: initialize state variables
x_dce_dDe = zeros(2,1); % initialize state variables
x_dce_depsilone = zeros(2,1); % initialize state variables
x_dce_depsilonesep = zeros(2,1); % initialize state variables
i_end = length(time); % final iteration
% Run SPMe
for i = 1:i_end
    [x_cse_n_updated_step, x_cse_p_updated_step, x_ce_updated_step, x_dcse_dDs_n_updated, x_dcse_dDs_p_updated, x_dcse_depsilons_n_updated, x_dcse_depsilons_p_updated, x_dce_dDe_updated, x_dce_depsilone_updated, x_dce_depsilonesep_updated, V_out_step(i), dVout_dDs_n_step(i), dVout_dDs_p_step(i), dVout_depsilons_n_step(i), dVout_depsilons_p_step(i), dVout_dk_n_step(i), dVout_dk_p_step(i), dVout_dDe_step(i), dVout_depsilone_step(i), dVout_depsilonesep_step(i), breakFlag] = SPMeSens_DS(I(i),D_s_p_true,epsilon_s_p_true,x_cse_n_step,x_cse_p_step,x_ce_step,x_dcse_dDs_n,x_dcse_dDs_p,x_dcse_depsilons_n,x_dcse_depsilons_p,x_dce_dDe,x_dce_depsilone,x_dce_depsilonesep,dt,V_out_max,V_out_min,SOC_max,SOC_min);
    x_cse_n_step = x_cse_n_updated_step; % update states
    x_cse_p_step = x_cse_p_updated_step;
    x_ce_step = x_ce_updated_step;
    x_dcse_dDs_n = x_dcse_dDs_n_updated;
    x_dcse_dDs_p = x_dcse_dDs_p_updated;
    x_dcse_depsilons_n = x_dcse_depsilons_n_updated;
    x_dcse_depsilons_p = x_dcse_depsilons_p_updated;
    x_dce_dDe = x_dce_dDe_updated;
    x_dce_depsilone = x_dce_depsilone_updated;
    x_dce_depsilonesep = x_dce_depsilonesep_updated;
    if breakFlag % stop condition
        i_end = i; % store final iteration (for indexing time vector when plotting)
        break
    end
end
% Store Values
true.D_s_p = D_s_p_true;
true.I = I;
true.SOC_0 = SOC_0;
true.V_out = V_out_step;
true.epsilon_s_p = epsilon_s_p_true;
true.norm_dVout_dDs_p = D_s_p_true*dVout_dDs_p_step;
true.norm_dVout_depsilons_p = epsilon_s_p_true*dVout_depsilons_p_step;
true.time = time;

%%% Run SPMe under Uncertain (a priori) Parameters
% Apply Uncertainty to Parameters
D_s_p_apriori = D_s_p_true*(1+uncertainty_Dsp/100); % [m^2/s] solid-phase diffusion coefficient for positive electrode
epsilon_s_p_apriori = epsilon_s_p_true*(1+uncertainty_epsilonsp/100); % [-] active material volume fraction of solid positive electrode
a_s_n = 3*epsilon_s_n/R_s_n; % [1/m] negative electrode specific interfacial surface area
a_s_p = 3*epsilon_s_p_apriori/R_s_p; % [1/m] positive electrode specific interfacial surface area
% Initialize
x_cse_n_step = [cse_0_n/(-10395*D_s_n^2/(a_s_n*F*R_s_n^5)); 0; 0]; % negative electrode: initial condition for state variables (uniform concentration profile)
x_cse_p_step = [cse_0_p/(-10395*D_s_p_apriori^2/(a_s_p*F*R_s_p^5)); 0; 0]; % positive electrode: initial condition for state variables (uniform concentration profile)
x_ce_step = 0; % state variable initial condition
x_dcse_dDs_n = zeros(4,1); % negative electrode: initialize state variables
x_dcse_dDs_p = zeros(4,1); % positive electrode: initialize state variables
x_dcse_depsilons_n = zeros(3,1); % negative electrode: initialize state variables
x_dcse_depsilons_p = zeros(3,1); % positive electrode: initialize state variables
x_dce_dDe = zeros(2,1); % initialize state variables
x_dce_depsilone = zeros(2,1); % initialize state variables
x_dce_depsilonesep = zeros(2,1); % initialize state variables
i_end = length(time); % final iteration
% Run SPMe
for i = 1:i_end
    [x_cse_n_updated_step, x_cse_p_updated_step, x_ce_updated_step, x_dcse_dDs_n_updated, x_dcse_dDs_p_updated, x_dcse_depsilons_n_updated, x_dcse_depsilons_p_updated, x_dce_dDe_updated, x_dce_depsilone_updated, x_dce_depsilonesep_updated, V_out_step(i), dVout_dDs_n_step(i), dVout_dDs_p_step(i), dVout_depsilons_n_step(i), dVout_depsilons_p_step(i), dVout_dk_n_step(i), dVout_dk_p_step(i), dVout_dDe_step(i), dVout_depsilone_step(i), dVout_depsilonesep_step(i), breakFlag] = SPMeSens_DS(I(i),D_s_p_apriori,epsilon_s_p_apriori,x_cse_n_step,x_cse_p_step,x_ce_step,x_dcse_dDs_n,x_dcse_dDs_p,x_dcse_depsilons_n,x_dcse_depsilons_p,x_dce_dDe,x_dce_depsilone,x_dce_depsilonesep,dt,V_out_max,V_out_min,SOC_max,SOC_min);
    x_cse_n_step = x_cse_n_updated_step; % update states
    x_cse_p_step = x_cse_p_updated_step;
    x_ce_step = x_ce_updated_step;
    x_dcse_dDs_n = x_dcse_dDs_n_updated;
    x_dcse_dDs_p = x_dcse_dDs_p_updated;
    x_dcse_depsilons_n = x_dcse_depsilons_n_updated;
    x_dcse_depsilons_p = x_dcse_depsilons_p_updated;
    x_dce_dDe = x_dce_dDe_updated;
    x_dce_depsilone = x_dce_depsilone_updated;
    x_dce_depsilonesep = x_dce_depsilonesep_updated;
    if breakFlag % stop condition
        i_end = i; % store final iteration (for indexing time vector when plotting)
        break
    end
end
% Store Values
apriori.D_s_p = D_s_p_apriori;
apriori.I = I;
apriori.SOC_0 = SOC_0;
apriori.V_out = V_out_step;
apriori.epsilon_s_p = epsilon_s_p_apriori;
apriori.norm_dVout_dDs_p = D_s_p_apriori*dVout_dDs_p_step;
apriori.norm_dVout_depsilons_p = epsilon_s_p_apriori*dVout_depsilons_p_step;
apriori.time = time;

%% Compute A Priori Rating Values
% A priori rating values (Q) are computed for each data segment in the full
% data set and organized in the 'Q-table'

% Extract Sensitivities
norm_dVout_dDs_p = apriori.norm_dVout_dDs_p; % [V] extract normalized D_s_p sensitivity
norm_dVout_depsilons_p = apriori.norm_dVout_depsilons_p; % [V] extract normalized epsilon_s_p sensitivity

% Load or Build Q-Table
if loadQtable % if loading Q-Table
    load(profile_filename+"_"+"Qtable_Dsp_Dsp"+uncertainty_Dsp+"_epsilonsp"+uncertainty_epsilonsp+".mat"); % load Q-Table
else % compute Q-Table
    tic
    Qtable = inf*ones(length(time)); % initialize Qtable (as table of infinies because we are interested in minimum values)
    loading = waitbar(0,"Creating Q Table");
    for i = 1:length(time) % for each data window size
        window = 2-i:1; % define data window (vector of indices)
        for ii = 1:length(time) % for each time step
            if window(1) >= 1 % if entire window is in data set
                sumSquareSens = norm_dVout_dDs_p(window)*norm_dVout_dDs_p(window)'; % compute denominator term of Q (sum of squared sensitivities)
                if sumSquareSens ~= 0 % if denominator of Q is not zero
                    Qtable(i,ii) = abs(norm_dVout_dDs_p(window)*norm_dVout_depsilons_p(window)')/sumSquareSens; % [-] compute data qualification rating
                end
            end
            window = window + 1; % increment window for next time step
        end
        waitbar(i/length(time));
    end
    close(loading);
    toc
    save(profile_filename+"_"+"Qtable_Dsp_Dsp"+uncertainty_Dsp+"_epsilonsp"+uncertainty_epsilonsp+".mat",'Qtable')
end
uniqueQ = unique(Qtable); % create array of ascending unique Q-values

%% Estimation Setup
% Setup
Voffset = 0; % [V] voltage offset to emulate measurement bias
V_true = true.V_out; % [V] extract voltage data
numTimeSteps = length(time); % number of time steps

% Estimation Parameters
paramUpperBound = 3e-13; % upper bound of estimated parameter
paramLowerBound = 1e-15; % lower bound of estimated parameter

% Apply Voltage Offset:
V_true = V_true + Voffset; % [V] apply voltage offset
        
% Load or Build Error Table
if loadDspErrorTable % if loading Dsp Error Table
    load(profile_filename+"_"+"DspVector_Dsp"+uncertainty_Dsp+"_epsilonsp"+uncertainty_epsilonsp+"_"+numDsp+".mat");
    load(profile_filename+"_"+"DspErrorTable_Dsp"+uncertainty_Dsp+"_epsilonsp"+uncertainty_epsilonsp+"_"+numDsp+".mat");
else % build Dsp Error Table
    DspVector = paramLowerBound:((paramUpperBound-paramLowerBound)/(numDsp-1)):paramUpperBound; % [m^2/s] Dsp values to run in SPMe
    DspVtable = zeros(length(DspVector),numTimeSteps); % initialize table for model voltage values under varying Dsp
    DspErrorTable = zeros(length(DspVector),numTimeSteps); % initialize table for squared voltage errors under varying Dsp
    tic
    loading = waitbar(0,"Creating Error Table");
    for ii = 1:length(DspVector) % for each value of Dsp
        a_s_n = 3*epsilon_s_n/R_s_n; % [1/m] negative electrode specific interfacial surface area
        a_s_p = 3*epsilon_s_p_apriori/R_s_p; % [1/m] positive electrode specific interfacial surface area        
        % Initialize States:
        x_cse_n_step = [cse_0_n/(-10395*D_s_n^2/(a_s_n*F*R_s_n^5)); 0; 0]; % initialize state variables (initial conditions from SPMe model)
        x_cse_p_step = [cse_0_p/(-10395*DspVector(ii)^2/(a_s_p*F*R_s_p^5)); 0; 0]; % initialize state variables with estimated parameter (initial conditions from SPMe model)
        x_ce_step = 0; % initialize state variable (initial condition from SPMe model)
        for i = 1:numTimeSteps
            [x_cse_n_updated_step, x_cse_p_updated_step, x_ce_updated_step, V_model] = SPMe(I(i),x_cse_n_step,x_cse_p_step,x_ce_step,DspVector(ii),epsilon_s_p_apriori,dt); % compute single time step of SPMe model
            x_cse_n_step = x_cse_n_updated_step; % append column of state variables
            x_cse_p_step = x_cse_p_updated_step; % append column of state variables
            x_ce_step = x_ce_updated_step; % append column of state variables
            DspVtable(ii,i) = V_model; % [V] save model voltage
            DspErrorTable(ii,i) = (V_model-V_true(i))^2; % [V^2] compute squared voltage error 
        end
        waitbar(ii/length(DspVector));
    end
    close(loading);
    toc
    save(profile_filename+"_"+"DspVector_Dsp"+uncertainty_Dsp+"_epsilonsp"+uncertainty_epsilonsp+"_"+numDsp+".mat",'DspVector');
    save(profile_filename+"_"+"DspErrorTable_Dsp"+uncertainty_Dsp+"_epsilonsp"+uncertainty_epsilonsp+"_"+numDsp+".mat",'DspErrorTable');
end

%% Perform Estimation and Compute A Posteriori Rating
% Sweep through data segments according to their a priori rating values to
% perform the estimation and compute the a posteriori rating value.

Qindex = 1:500; % [-] indices of Q-table associated with data segments to be used for estimation (i.e., Qindex = 1:500 will sweep through the 500 data segments with the smallest a priori rating values to perform the estimation and compute the a posteriori rating value)
counter = 1; % initialize counter

for iteration = 1:length(Qindex)
    [window_length,window_lastTimeStep] = find(Qtable == uniqueQ(Qindex(iteration))); % look up length and final time step of the selected data segment
    selectedData = [window_lastTimeStep-window_length+1:window_lastTimeStep]; % indices of selected data
    Q_vector(counter) = uniqueQ(Qindex(iteration)); % [-] store a priori rating value
    
    % Perform Estimation
    sumSquareError = sum(DspErrorTable(:,selectedData),2)'; % [V^2] compute sum of squared error
    Dsp_estimate = DspVector(sumSquareError==min(sumSquareError)); % [m^2/s] Dsp estimate (from tabular sum of squares estimation)
    
    %%% Run SPMe under Estimated D_s_p
    a_s_n = 3*epsilon_s_n/R_s_n; % [1/m] negative electrode specific interfacial surface area
    a_s_p = 3*epsilon_s_p_apriori/R_s_p; % [1/m] positive electrode specific interfacial surface area
    % Initialize
    x_cse_n_step = [cse_0_n/(-10395*D_s_n^2/(a_s_n*F*R_s_n^5)); 0; 0]; % negative electrode: initial condition for state variables (uniform concentration profile)
    x_cse_p_step = [cse_0_p/(-10395*Dsp_estimate^2/(a_s_p*F*R_s_p^5)); 0; 0]; % positive electrode: initial condition for state variables (uniform concentration profile)
    x_ce_step = 0; % state variable initial condition
    x_dcse_dDs_n = zeros(4,1); % negative electrode: initialize state variables
    x_dcse_dDs_p = zeros(4,1); % positive electrode: initialize state variables
    x_dcse_depsilons_n = zeros(3,1); % negative electrode: initialize state variables
    x_dcse_depsilons_p = zeros(3,1); % positive electrode: initialize state variables
    x_dce_dDe = zeros(2,1); % initialize state variables
    x_dce_depsilone = zeros(2,1); % initialize state variables
    x_dce_depsilonesep = zeros(2,1); % initialize state variables
    i_end = length(time); % final iteration
    % Run SPMe
    for i = 1:i_end
        [x_cse_n_updated_step, x_cse_p_updated_step, x_ce_updated_step, x_dcse_dDs_n_updated, x_dcse_dDs_p_updated, x_dcse_depsilons_n_updated, x_dcse_depsilons_p_updated, x_dce_dDe_updated, x_dce_depsilone_updated, x_dce_depsilonesep_updated, V_out_step(i), dVout_dDs_n_step(i), dVout_dDs_p_step(i), dVout_depsilons_n_step(i), dVout_depsilons_p_step(i), dVout_dk_n_step(i), dVout_dk_p_step(i), dVout_dDe_step(i), dVout_depsilone_step(i), dVout_depsilonesep_step(i), breakFlag] = SPMeSens_DS(I(i),Dsp_estimate,epsilon_s_p_apriori,x_cse_n_step,x_cse_p_step,x_ce_step,x_dcse_dDs_n,x_dcse_dDs_p,x_dcse_depsilons_n,x_dcse_depsilons_p,x_dce_dDe,x_dce_depsilone,x_dce_depsilonesep,dt,V_out_max,V_out_min,SOC_max,SOC_min);
        x_cse_n_step = x_cse_n_updated_step; % update states
        x_cse_p_step = x_cse_p_updated_step;
        x_ce_step = x_ce_updated_step;
        x_dcse_dDs_n = x_dcse_dDs_n_updated;
        x_dcse_dDs_p = x_dcse_dDs_p_updated;
        x_dcse_depsilons_n = x_dcse_depsilons_n_updated;
        x_dcse_depsilons_p = x_dcse_depsilons_p_updated;
        x_dce_dDe = x_dce_dDe_updated;
        x_dce_depsilone = x_dce_depsilone_updated;
        x_dce_depsilonesep = x_dce_depsilonesep_updated;
        if breakFlag % stop condition
            i_end = i; % store final iteration (for indexing time vector when plotting)
            break
        end
    end
    % Compute A Posteriori Quality Rating
    actualError = (Dsp_estimate-true.D_s_p)/true.D_s_p; % [-] true estimation error
    Q_underEstimatedParameters = abs((Dsp_estimate*dVout_dDs_p_step(selectedData))*(apriori.epsilon_s_p*dVout_depsilons_p_step(selectedData)'))/sum((Dsp_estimate*dVout_dDs_p_step(selectedData)).^2); % [-] compute a posteriori data quality rating

    % Store Values
    AE_vector(counter) = actualError;
    BF_est(counter) = Dsp_estimate;
    windowLength_vector(counter) = length(selectedData);
    windowLastTimeStep_vector(counter) = selectedData(end);
    Q_underEstimatedParameters_vector(counter) = Q_underEstimatedParameters;
    
    counter = counter+1; % increment counter
end

%% Final Estimation
estGroupThreshold = 0.05; % [-] fraction of beginning values to use for the estimation
estMatrix = [Q_underEstimatedParameters_vector' BF_est']; % build matrix for estimation: column 1 = a posteriori ratings under estimated parameters; column 2 = estimates
estMatrix = sortrows(estMatrix,1); % sort rows by ascending a posteriori ratings
estGroupThresholdIndex = round(estGroupThreshold*length(BF_est)); % compute cut-off index for the estimation group
Dsp_estimate_final_avg = mean(estMatrix(1:estGroupThresholdIndex,2)) % mean parameter estimate from the subset of segments with the smallest a posteriori ratings
sumSquareError_fullprofile = sum(DspErrorTable(:,:),2)'; % [V^2] compute sum of squared error for full profile
Dsp_estimate_fullprofile = DspVector(sumSquareError_fullprofile==min(sumSquareError_fullprofile)) % [m^2/s] Dsp estimate under full profile (no data selection)

%% Output Results
results.Qindex = Qindex';
results.Qapriori = Q_vector';
results.estimate = BF_est';
results.Qaposteriori = Q_underEstimatedParameters_vector';
results.ActualError = AE_vector';
results = struct2table(results);

%% Plots
fontSize = 11;
plot_index1 = 20; % a priori quality rating index to plot as an example segment
figure('Renderer', 'painters', 'Position', [680,558-250,560,420*1.5]);

% Current Profile
subplot(3,1,1)
selectDataSeg1 = [windowLastTimeStep_vector(plot_index1)-windowLength_vector(plot_index1)+1:windowLastTimeStep_vector(plot_index1)]; % indices of selected data
plot(time,I,'k',time(selectDataSeg1),I(selectDataSeg1),'r')
xlabel('Time $(s)$','interpreter','latex','fontsize',fontSize);
ylabel('Discharge Current $(A)$','interpreter','latex','fontsize',fontSize);
title("Input Current Profile",'interpreter','latex','fontsize',fontSize)
legend("Input Data Set","Example Data Segment","Location","Best",'interpreter','latex','fontsize',fontSize-1)
set(gca,'TickLabelInterpreter','latex')
grid on

% Estimation Error vs. A Priori Quality Rating
subplot(3,1,2); hold on
scatter(Q_vector,abs(AE_vector)*100)
scatter(Q_vector(plot_index1),abs(AE_vector(plot_index1))*100,'filled','r')
xlabel('\textit{a priori} Data Quality Rating, $Q^{-}_{D_{s,p}}\ (-)$','interpreter','latex','fontsize',fontSize);
ylabel('Estimation Error $(\%)$','interpreter','latex','fontsize',fontSize);
title("$D_{s,p}$ Estimation Error vs. \textit{a priori} Data Quality Rating",'interpreter','latex','fontsize',fontSize)
legend("Data Segments","Example Data Segment","Location","Best",'interpreter','latex','fontsize',fontSize-1)
xlim([0 0.12])
set(gca,'TickLabelInterpreter','latex')
grid on

% Estimation Error vs. A Posteriori Quality Rating
subplot(3,1,3); hold on
scatter(Q_underEstimatedParameters_vector,abs(AE_vector)*100)
scatter(Q_underEstimatedParameters_vector(plot_index1),abs(AE_vector(plot_index1))*100,'filled','r')
xlabel('\textit{a posteriori} Data Quality Rating, $Q^{+}_{D_{s,p}}\ (-)$','interpreter','latex','fontsize',fontSize);
ylabel('Estimation Error $(\%)$','interpreter','latex','fontsize',fontSize);
title("$D_{s,p}$ Estimation Error vs. \textit{a posteriori} Data Quality Rating",'interpreter','latex','fontsize',fontSize)
legend("Data Segments","Example Data Segment","Location","Best",'interpreter','latex','fontsize',fontSize-1)
xlim([0 6])
set(gca,'TickLabelInterpreter','latex')
grid on
