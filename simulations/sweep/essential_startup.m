%% Declare and add essential paths
clc; clear; close all;
head_dir = '../..';

prm_folder  = fullfile(head_dir, 'prm');
func_folder = fullfile(head_dir, 'functions');

ext_pck2add = {'yamlmatlab', 'GetFullPath'}; 
ext_folder = fullfile(head_dir, 'extpckgs'); 

addpath(genpath(prm_folder));
addpath(genpath(func_folder));

cellfun(@(x) addpath(genpath(fullfile(ext_folder, x))), ext_pck2add); 
