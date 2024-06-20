# Active_RL_for_FI_objective
# Rui Huang
# 2024-06-16

This is the repo for actitive reinforcement learning code with Fisher information objective to generate current excitation for estimating EPS battery parameters. Refering to paper "Reinforcement learning of optimal input excitation for parameter estimation with application to li-ion battery.", Huang, Rui, Jackson Fogelquist, and Xinfan Lin. IEEE Transactions on Industrial Informatics (2023).

It includes:

(1) QLearning_main.m: the main file serving as an example to train the RL agent and save the excitation optimization results

(2) SPMe_Step_fcn.m: function to propagate the battery model and serve as the RL environment

(3) parameters.ini: define the values of battery parameters  

(4) inifile.m: creates, reads, or writes data from/to a standard ini (ascii) file.

(5) IniConfig.m The class for working with configurations of settings and INI-files.

Input:
No inputs required

Output:
Current profile time-series data, saved in the "learned_policy.csv" in the root directory.
Reinforcement learning policy, saved in the "Q_table.csv" in the root directory.

Code written and conducted on a Matlab R2021a platfrom on a Windows 10 machine

