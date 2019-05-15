
%% Load set up files and add path 
load set_up.mat; 
load sim_var.mat; 

%% Define "edge_lim" for rate calculation 
set_up.rate_estimation.edge_lim = [0, set_up.time.tstop]*1e-3; 

%% Create variable combination and *SIMGLOB 
connectivity_weight_EtoE_k_SIMGLOB = sim_var.connectivity_weight_EtoE_k; 
connectivity_conductance_EtoE_SIMGLOB = sim_var.connectivity_conductance_EtoE; 
connectivity_conductance_EtoI_SIMGLOB = sim_var.connectivity_conductance_EtoI; 
[combo_vec_SIMGLOB, ~, combo_struct_SIMGLOB] = return_combination( ... 
	connectivity_weight_EtoE_k_SIMGLOB, ...
	connectivity_conductance_EtoE_SIMGLOB, ...
	connectivity_conductance_EtoI_SIMGLOB); 
num_SIMGLOB = size(combo_vec_SIMGLOB,1); 
combo_struct_SIMGLOB.files = cell(num_SIMGLOB, 1); 

%% Sweep loop 
for i_SIMGLOB = 1 : num_SIMGLOB 
	 file_prefix = sprintf('sweep_%03d', i_SIMGLOB); 
	 combo_struct_SIMGLOB.files{i_SIMGLOB} = file_prefix; 
	 vec_simvar = combo_vec_SIMGLOB(i_SIMGLOB,:); 

	%% Change variables within `set_up`
	 set_up.file_prefix = file_prefix; 

	% Define "num" for stimulation purposes
	 set_up.stimulation.num = ceil(set_up.stimulation.portion*set_up.neurons.num); 

	%% Change requested variations
	 set_up.connectivity.weight.EtoE.k = vec_simvar{1}; 
	 set_up.connectivity.conductance.EtoE = vec_simvar{2}; 
	 set_up.connectivity.conductance.EtoI = vec_simvar{3}; 

	%% INITIALIZE SIMULATION
	 clearvars -except set_up *SIMGLOB; 
	 run run_sim.m; 
	 close all; 
end 

%% Save `sweep_master` for book-keeping purposes
sweep_master = combo_struct_SIMGLOB; 
save(fullfile(set_up.paths.data, 'sweep_master.mat'), 'sweep_master'); 
