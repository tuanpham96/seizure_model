run essential_startup.m

import setup_sweep.*
yaml_setup_file = 'set_up.yaml'; 
sim_path = 'sweepE2E'; 
tmpl_path = 'template_scripts'; 
rel_sim_data_path = '..'; 
tic
parse_setupfile(yaml_setup_file, sim_path, tmpl_path, rel_sim_data_path); 
toc

run(fullfile(sim_path, 'init_sweep.m'))