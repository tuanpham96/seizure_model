function parse_setupfile(yaml_setup_file, sim_path)
if nargin == 1
    sim_path = '.'; 
end

import yaml.*
struct_obj = ReadYaml(yaml_setup_file); 
field_list = recursive_fieldnames(struct_obj); 
[field_names, var_vecs, var_labels, var_alias] = find_sim(struct_obj, field_list); 

glob_pref = return_field_value(struct_obj, 'alias', 'SWEEP'); 

copyfile(yaml_setup_file, fullfile(sim_path, yaml_setup_file)); 
init_script = fopen(fullfile(sim_path, 'init_sweep.m', 'w');
set_script = fopen(fullfile(sim_path, 'set_sweep.m', 'w');

fprintf(init_script, '
end

function [field_names, var_vecs, var_labels, var_alias] = find_sim(struct_obj, field_list)
ind_with_vec = cellfun(@isempty, regexp(field_list, '.vec$', 'match'), 'uni', 1); 
vec_fields = field_list(ind_with_vec); 
var_vecs = cellfun(@(x) getfield_recurse(struct_obj, x), vec_fields, 'uni', 0);

field_names = celfun(@(x) x(1:end-4), vec_fields, 'uni', 0); 

lbl_fields = cellfun(@(x) [x, '.label'], field_names, 'uni', 0);
has_lbl = find(cellfun(@(x) strcmp(field_list, x), lbl_fields, 'uni', 1)); 

var_labels = field_names; 
var_labels(has_lbl) = cellfun(@(x) getfield_recurse(struct_obj, lbl_fields(x)), has_lbl, 'uni', 0); 

var_alias = strcat(regexprep(field_names, '.', '_'), '_SIMGLOB'); 

end