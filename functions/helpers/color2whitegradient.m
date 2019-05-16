function cmap = color2whitegradient(c, n)
cmap = arrayfun(@(ci) to_col_vec(linspace(ci, 1, n)), c, 'uni', 0); 
cmap = horzcat(cmap{:}); 
end