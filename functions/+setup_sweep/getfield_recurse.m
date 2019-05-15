function v = getfield_recurse(s, f)
f = strsplit(f, '.'); 
v = getfield(s, f{:});  
end