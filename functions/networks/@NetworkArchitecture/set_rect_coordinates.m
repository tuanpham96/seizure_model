function set_rect_coordinates(obj, rect_size)

% First obtain a grid
N = obj.num_neurons;
lenX = rect_size(1);
lenY = rect_size(2);

xvec = 1:lenX;
yvec = 1:lenY;

combo_coord = return_combomat(xvec, yvec);

% Shuffle the order
shuff_coord = combo_coord(randperm(N), :);

x = shuff_coord(:,1);
y = shuff_coord(:,2);

% Euclidean distance
dist_mat = euclidean_distance_matrix(x, y);

% Get center coordinate
xcent = ceil(lenX/2);
ycent = ceil(lenY/2);
cent_ind = find(x == xcent & y == ycent, 1);

% Get side coordinate 
xside = min(x); 
yside = ycent; 
side_ind = find(x == xside & y == yside, 1); 

% Get corner coordinate
xcorn = xside; 
ycorn = min(y); 
corn_ind = find(x == xcorn & y == ycorn, 1); 

% Save to struct and fields
obj.coord = struct( 'x', x, 'y', y, 'size', rect_size);
obj.dist = dist_mat;
obj.centered_neuron = struct( 'x', xcent, 'y', ycent, ...
    'ind', cent_ind);
obj.side_neuron = struct( 'x', xside, 'y', yside, 'ind', side_ind); 
obj.corner_neuron = struct( 'x', xcorn, 'y', ycorn, 'ind', corn_ind); 
end