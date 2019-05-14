clc; clear; close all; 
addpath(genpath('../../functions')); 
addpath(genpath('../../extpckgs/yamlmatlab')); 

import setup_sweep.*
yaml_setup_file = 'test.yaml'; 
sim_path = 'sweep_test'; 

tic
parse_setupfile(yaml_setup_file, sim_path)
toc