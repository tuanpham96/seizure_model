%% Set Up
set_up.file_pref = 'Test 1';
%% Specific Matrix
set_up.FFW = [0,1,1;0,0,0;0,-1,0]; 
set_up.FB = [0,1,0;0,0,1;0,-1,0];
%% Test Matrix
set_up.TM = zeros(100,100);
%% Random Subgraph Sampling
set_up.percent = 0.02;
%% Iterations
set_up.iter = 10;