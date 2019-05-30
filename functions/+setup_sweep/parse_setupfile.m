function parse_setupfile(yaml_setup_file, sim_path, sim_tmpl_path, rel_sim_data_path)

import yaml.* %#ok<NSTIMP>
import setup_sweep.*  %#ok<NSTIMP>

% Read YAML file and obtain:
%   set_up:         cleaned `set_up` struct
%   sim_var:        variable structs (meaning the ones being varied) 
%   field_names:    the actual field levels to get variables 
%   var_alias:      the aliases to name the variable in script 
struct_obj = ReadYaml(yaml_setup_file);
field_list = recursive_fieldnames(struct_obj);
[set_up, sim_var, field_names, var_alias] = find_sim_var(struct_obj, field_list);

glob_pref = return_field_value(struct_obj, 'file_prefix', 'sweep');
if ~exist(sim_path, 'dir') 
    mkdir(sim_path);
end

    function path_name = sim_file_path(file_name)
        path_name = fullfile(sim_path, file_name);
    end

% Make "raw" and "processed" data folders + script folder
% the script folder is just for referencing - currently the scripts to run
% the simulation is stored in `sim_path` for version control purposes;
% whereas the script folder in the data folder is just another copy for
% referencing 

set_up.paths.data = fullfile(set_up.paths.data, set_up.file_prefix); 
paths = struct(); 
paths.data = set_up.paths.data; 
paths.raw_data = fullfile(paths.data, 'raw'); 
paths.processed_data = fullfile(paths.data, 'processed');
paths.script = fullfile(paths.data, 'script'); 

path_fields = fieldnames(paths); 
for i = 1:length(path_fields) 
    field_i = path_fields{i}; 
    path_i = paths.(field_i); 
    
    % create folder if not already there 
    if ~exist(path_i, 'dir') 
        mkdir(path_i);
    end
    
    % this is for referencing when saving 
    set_up.paths.(field_i) = fullfile(rel_sim_data_path, path_i);
end

% Save structs, yaml files (book-keeping) and initialization script 
save(sim_file_path('set_up.mat'), 'set_up'); 
save(sim_file_path('sim_var.mat'), 'sim_var'); 
copyfile(yaml_setup_file, sim_file_path('set_up.yaml'));
init_script = fopen(sim_file_path('init_sweep.m'), 'w');

% Copy all templates script to folder
copyfile([sim_tmpl_path '/*'], sim_path); 

% Load set up files 
fprintf(init_script, '\n%%%% Load set up files and add path \n'); 
fprintf(init_script, 'load set_up.mat; \n'); 
fprintf(init_script, 'load sim_var.mat; \n'); 

% Define *SIMGLOB variables and create COMBO
fprintf(init_script, '\n%%%% Create variable combination and *SIMGLOB \n');
cellfun(@(alias) fprintf(init_script, '%s_SIMGLOB = sim_var.%s; \n', alias, alias), var_alias); 

sep_all_struct = sprintf(', ...\n'); 
all_struct_simglob = sprintf(['\t%s_SIMGLOB' sep_all_struct], var_alias{:});  
all_struct_simglob(end-length(sep_all_struct)+1:end) = ''; 

fprintf(init_script, ...
    '[combo_vec_SIMGLOB, ~, combo_struct_SIMGLOB] = return_combination( ... \n%s); \n', ...
    all_struct_simglob);

fprintf(init_script, [...
    'num_SIMGLOB = size(combo_vec_SIMGLOB,1); \n' ...
    'combo_struct_SIMGLOB.files = cell(num_SIMGLOB, 1); \n']);

% Write loop to simulate 
fprintf(init_script, [...
    '\n%%%% Sweep loop \n' ...
    'for i_SIMGLOB = 1 : num_SIMGLOB \n' ...
    '\t file_prefix = sprintf(''%s_%%03d'', i_SIMGLOB); \n' ...
    '\t combo_struct_SIMGLOB.files{i_SIMGLOB} = file_prefix; \n' ...
    '\t vec_simvar = combo_vec_SIMGLOB(i_SIMGLOB,:); \n' ...
    '\n\t%%%% Change variables within `set_up`\n' ...
    ], glob_pref); 

fprintf(init_script, '\t set_up.file_prefix = file_prefix; \n'); 

% Changing variables 
fprintf(init_script, '\n\t%%%% Change requested variations\n'); 
arrayfun(@(i) fprintf(init_script, '\t set_up.%s = vec_simvar{%d}; \n', ...
    field_names{i}, i) , 1:length(field_names));

% Initialize simulation 
fprintf(init_script, [ ...
    '\n\t%%%% INITIALIZE SIMULATION\n' ...
    '\t clearvars -except set_up *SIMGLOB; \n' ...
    '\t run run_sim.m; \n' ...
    '\t close all; \n' ...
    'end \n']);

% Save sweep master 
fprintf(init_script, [ ...
    '\n%%%% Save `sweep_master` for book-keeping purposes\n' ...
    'sweep_master = combo_struct_SIMGLOB; \n' ... 
    'save(fullfile(set_up.paths.data, ''sweep_master.mat''), ''sweep_master''); \n']);  

% Copy simulation script folder to data folder 
copyfile([sim_path '/*'], paths.script); 

end

