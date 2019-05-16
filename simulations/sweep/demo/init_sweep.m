
%% Load set up files and add path 
load set_up.mat; 
load sim_var.mat; 

%% Create variable combination and *SIMGLOB 
stimulation_location_SIMGLOB = sim_var.stimulation_location; 
[combo_vec_SIMGLOB, ~, combo_struct_SIMGLOB] = return_combination( ... 
	stimulation_location_SIMGLOB); 
num_SIMGLOB = size(combo_vec_SIMGLOB,1); 
combo_struct_SIMGLOB.files = cell(num_SIMGLOB, 1); 

%% Sweep loop 
for i_SIMGLOB = 1 : num_SIMGLOB 
	 file_prefix = sprintf('sweep_%03d', i_SIMGLOB); 
	 combo_struct_SIMGLOB.files{i_SIMGLOB} = file_prefix; 
	 vec_simvar = combo_vec_SIMGLOB(i_SIMGLOB,:); 

	%% Change variables within `set_up`
	 set_up.file_prefix = file_prefix; 

	%% Change requested variations
	 set_up.stimulation.location = vec_simvar{1}; 

	%% INITIALIZE SIMULATION
	 clearvars -except set_up *SIMGLOB; 
	 run run_sim.m; 
	 close all; 
end 

%% Save `sweep_master` for book-keeping purposes
sweep_master = combo_struct_SIMGLOB; 
save(fullfile(set_up.paths.data, 'sweep_master.mat'), 'sweep_master'); 
