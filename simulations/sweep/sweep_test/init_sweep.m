
%% Load set up files and add path 
load set_up.mat; 
load sim_var.mat; 
addpath(genpath('sweep_test')); 
addpath(genpath(set_up.paths.functions)); 
addpath(genpath(set_up.paths.parameters)); 
addpath(genpath(set_up.paths.data)); 

%% Create variable combination and *SIMGLOB 
weight_prm_EtoE_k_SIMGLOB = sim_var.weight_prm_EtoE_k; 
weight_prm_EtoI_type_SIMGLOB = sim_var.weight_prm_EtoI_type; 
weight_prm_EtoI_dist_range_SIMGLOB = sim_var.weight_prm_EtoI_dist_range; 
Gsyn_EtoE_SIMGLOB = sim_var.Gsyn_EtoE; 
Gsyn_EtoI_SIMGLOB = sim_var.Gsyn_EtoI; 
stim_percent_stim_SIMGLOB = sim_var.stim_percent_stim; 
rando_vec_lvl_1_SIMGLOB = sim_var.rando_vec_lvl_1; 
[combo_vec_SIMGLOB, ~, combo_struct_SIMGLOB] = return_combination( ... 
	weight_prm_EtoE_k_SIMGLOB, ...
	weight_prm_EtoI_type_SIMGLOB, ...
	weight_prm_EtoI_dist_range_SIMGLOB, ...
	Gsyn_EtoE_SIMGLOB, ...
	Gsyn_EtoI_SIMGLOB, ...
	stim_percent_stim_SIMGLOB, ...
	rando_vec_lvl_1_SIMGLOB); 
num_SIMGLOB = size(combo_vec_SIMGLOB,1); 
combo_struct_SIMGLOB.files = cell(num_SIMGLOB, 1); 

%% Sweep loop 
for i_SIMGLOB = 1 : num_SIMGLOB 
	 file_pref = sprintf('sweep_me_now_%03d', i_SIMGLOB); 
	 combo_struct_SIMGLOB.files{i_SIMGLOB} = file_pref; 
	 vec_simvar = combo_vec_SIMGLOB(i_SIMGLOB,:); 

	%% Change variables within `set_up`
	 set_up.file_pref = file_pref; 

	% Define "num_stimulated" for stimulation purposes
	 set_up.stim.num_stimulated = ceil(set_up.stim.percent_stim*set_up.num_neurons); 

	%% Change requested variations
	 set_up.weight_prm.EtoE.k = vec_simvar{1}; 
	 set_up.weight_prm.EtoI.type = vec_simvar{2}; 
	 set_up.weight_prm.EtoI.dist_range = vec_simvar{3}; 
	 set_up.Gsyn.EtoE = vec_simvar{4}; 
	 set_up.Gsyn.EtoI = vec_simvar{5}; 
	 set_up.stim.percent_stim = vec_simvar{6}; 
	 set_up.rando_vec.lvl_1 = vec_simvar{7}; 

	%% INITIALIZE SIMULATION
	 clearvars -except set_up *SIMGLOB; 
	 run run_sim.m; 
	 close all; 
end 

%% Save `sweep_master` for book-keeping purposes
sweep_master = combo_struct_SIMGLOB; 
save(fullfile(set_up.paths.data, 'sweep_master.mat'), 'sweep_master'); 
