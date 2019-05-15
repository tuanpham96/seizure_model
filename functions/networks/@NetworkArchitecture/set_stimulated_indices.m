function stim_ind = set_stimulated_indices(obj, location, num_stim) 
switch upper(location)
    case 'CENTER'
        prop_name = 'centered_neuron'; 
    case 'SIDE' 
        prop_name = 'side_neuron';
    case 'CORNER' 
        prop_name = 'corner_neuron'; 
    otherwise 
        error('The input "location" = %s is not allowed. Only CENTER, SIDE or CORNER', location);
end

stim_nrn = obj.(prop_name); 
if isempty(stim_nrn)
    error(['Cannot proceed, possibly because the spatial distribution and ' ...
        'the stimulation location are incompatible']); 
end

start_ind = stim_nrn.ind; 
neigh_ind = obj.pick_random_neighbors(num_stim-1,start_ind); 
stim_ind = [start_ind, neigh_ind]; 


end